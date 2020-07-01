#!/usr/bin/env python
# coding=utf-8
"""
transcript_expression.py

Part of neoepiscope
Functions for determining expression of transcripts.

Licensed under the MIT license.

The MIT License (MIT)
Copyright (c) 2018 Mary A. Wood, Austin Nguyen,
                   Abhinav Nellore, and Reid Thompson

Portions from Rail-RNA, which is copyright (c) 2015 
                    Abhinav Nellore, Leonardo Collado-Torres,
                    Andrew Jaffe, James Morton, Jacob Pritt,
                    José Alquicira-Hernández,
                    Christopher Wilks,
                    Jeffrey T. Leek, and Ben Langmead.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import absolute_import, division, print_function
import re

def feature_to_tpm_dict(feature_to_read_count, feature_to_feature_length):
    """ Calculate TPM values for feature

        feature_to_read_count: dictionary linking features to read counts (float)
        feature_to_feature_length: dictionary linking features to feature lengths (float)

        Return value: dictionary linking feature ID to TPM value
    """
    total_rpk = 0.0
    feature_to_rpk = {}
    feature_to_tpm = {}
    # Get read per kilobase counts for each feature
    for feature in feature_to_read_count:
        try:
            rpk = feature_to_read_count[feature]/feature_to_feature_length[feature]
        except KeyError:
            continue
        feature_to_rpk[feature] = rpk
        total_rpk += rpk
    # Calculate scaling factor
    scaling = total_rpk/1000000.0
    # Calculate TPM values
    for feature in feature_to_rpk:
        tpm = feature_to_rpk[feature]/scaling
        feature_to_tpm[feature] = tpm
    return feature_to_tpm

# Taken from Rail-RNA: https://github.com/nellore/rail
def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers
        md: an MD string (example: 33A^CC).
        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

# Taken from Rail-RNA: https://github.com/nellore/rail
def indels_junctions_exons_mismatches(cigar, md, pos, seq, drop_deletions=False, junctions_only=False):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list
        Return value: tuple (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion, 
            string of inserted bases). Deletions is a list of tuples (first genomic 
            position of deletion, string of deleted bases). Junctions is a list
            of tuples (intron start position (inclusive), intron end position (exclusive),
            left_diplacement, right_displacement). Exons is a list of tuples (exon start 
            position (inclusive), exon end position (exclusive)). Mismatches is a list
            of tuples (genomic position of mismatch, read base)
    """
    insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(['cigar and md:', ''.join(cigar), ''.join(md)])
                    if not junctions_only:
                        mismatches.append((pos + aligned_bases, seq[seq_index + aligned_bases]))
                    correction_length = len(md[md_index])
                    m_length = aligned_base_cap - aligned_bases
                    if correction_length > m_length:
                        md[md_index] = md[md_index][:m_length]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += correction_length
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add junction
            junctions.append((pos, pos + skip_increment, seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append((pos - 1, seq[seq_index:seq_index+insert_size]))
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(['cigar and md:', ''.join(cigar), ''.join(md)])
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if not drop_deletions: exons.append((pos, pos + delete_size))
            if md_delete_size > delete_size:
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            assert cigar[cigar_index+1] == 'S'
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    # Merge exonic chunks/deletions; insertions/junctions could have chopped them up
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
        new_exons.append(last_exon)
    return insertions, deletions, junctions, new_exons, mismatches

