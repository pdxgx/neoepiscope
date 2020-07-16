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
from intervaltree import IntervalTree
from collections import defaultdict
import re
import os
import ast
import copy
import pysam
import tempfile
import subprocess

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
def indels_junctions_exons_mismatches(
            cigar, md, pos, seq, drop_deletions=False, junctions_only=False
        ):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS
    
        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list

        Return value: tuple
            (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Junctions is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)). Mismatches is a list
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
                    assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
                    if not junctions_only:
                        mismatches.append(
                                (pos + aligned_bases,
                                    seq[seq_index + aligned_bases])
                            )
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
            junctions.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if not drop_deletions: exons.append((pos, pos + delete_size))
            if md_delete_size > delete_size:
                # Deletion contains a junction
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index+1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/junctions could have chopped
    them up.'''
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

def generate_variant_bed(neopeptide_dict, chr_in_contigs):
    ''' Generates bed file with locations of neopeptide-causing variants

        neopeptide_dict: dictionary linking neopeptide sequences to list of
                         annotation information, where each annotation
                         consists of (chromosome, position, reference allele,
                         alternative allele, variant type, VAF, paired normal
                         epitope, transcript warnings, transcript identifier)
        chr_in_contigs: whether bam contigs contain 'chr' prefix (boolean)

        Return value: path to temporary bed file (str), 
                      set of neopeptide-causing variants,
                      dictionary linking chromosome to interval tree with
                      intervals linked to variants
    '''
    # Get set of all neopeptide-causing variants
    all_mutations = set()
    for pep in neopeptide_dict:
        mutations = set([(x[0], x[1], x[2], x[3], x[4]) for x in neopeptide_dict[pep]])
        all_mutations.update(mutations)
    # Create bed file for these variants
    bed_path = tempfile.mkstemp(suffix=".variant.bed", text=True)[1]
    var_intervals = defaultdict(IntervalTree)
    with open(bed_path, 'w') as f:
        for mut in all_mutations:
            # Figure out chromosome name structure
            contig = mut[0]
            if (mut[0].startswith('chr') and chr_in_contigs) or (not mut[0].startswith('chr') and not chr_in_contigs):
                contig = mut[0]
            elif mut[0].startswith('chr') and not chr_in_contigs:
                contig = mut[0].replace('chr', '')
            elif not mut[0].startswith('chr') and chr_in_contigs:
                contig = ''.join(['chr', mut[0]])
            # Find end positin of variant
            if mut[4] == 'V':
                ref = mut[2]
                alt = mut[3]
                end = mut[1] + len(mut[3])
            elif mut[4] == 'I':
                ref = ''
                alt = mut[3]
                end = mut[1] + len(mut[3]) + 1
            elif mut[4] == 'D':
                ref = mut[2]
                alt = mut[3]
                end = mut[1] + mut[3]
            out_line = [contig, str(mut[1]), str(end)]
            var_intervals[contig][mut[1]:end] = (mut[0], mut[1], ref, alt, mut[4])
            print('\t'.join(out_line), file=f)
    return bed_path, all_mutations, var_intervals

### Taken with modifications from https://www.biostars.org/p/306041/
def read_pair_generator(bam, contig):
    """ Generate read pairs in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.

        bam: path to BAM file
        contig: name of contig to process

        Return value: iterator of read1, read2 for each paired read
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(contig=contig):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            # Ignore reads that aren't mapped as pairs or from primary alignments
            continue
        if read.query[-2:] in ['.1', '.2']:
            # Strip '.1' or '.2' from query name
            qname = read.query_name[:-2]
        else:
            qname = read.query_name
        if qname not in read_dict:
            # Establishing new read pair
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            # Finish established read pair
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def get_expressed_variants(bam, reference_index, neopeptides, remove_files=True):
    ''' Gets transcripts that are expressed by at least one read pair

        bam: path to RNA-seq BAM file
        reference_index: bowtie reference index
        neopeptides: dictionary linking neopeptide sequences to list of
                     annotation information, where each annotation
                     consists of (chromosome, position, reference allele,
                     alternative allele, variant type, VAF, paired normal
                     epitope, transcript warnings, transcript identifier)

        Return value: set of transcript IDs that are expressed
    '''
    files_to_remove = []
    expressed_variants = defaultdict(int)
    covered_variants = defaultdict(int)
    # Check contig names
    chr_in_contigs = False
    orig_bam = pysam.AlignmentFile(bam, 'rb')
    if [x for x in orig_bam.references]:
        chr_in_contigs = True
    # Set temporary files
    new_bam = tempfile.mkstemp(suffix=".reduced.bam", text=True)[1]
    files_to_remove.append(new_bam)
    files_to_remove.append('.'.join([new_bam, 'bai']))
    bed_path, all_mutations, var_intervals = generate_variant_bed(neopeptides, chr_in_contigs)
    files_to_remove.append(bed_path)
    # Create and index reduced bam file
    subprocess.check_call(['samtools', 'view', '-b', '-L', bed_path, '-o', new_bam, bam])
    subprocess.check_call(['samtools', 'index', new_bam, ''.join([new_bam, '.bai'])])
    # Process BAM file
    bam_reader = pysam.AlignmentFile(new_bam, 'rb')
    for contig in reference_index.recs.keys():
        # ID contig name
        if contig in bam_reader.references:
            search_contig = copy.copy(contig)
        elif contig.replace('chr', '') in bam_reader.references:
            search_contig = contig.replace('chr', '')
        elif ''.join(['chr', contig]) in bam_reader.references:
            search_contig = ''.join(['chr', contig])
        else:
            continue
        for read1, read2 in read_pair_generator(bam_reader, search_contig):
            # Get read info
            r1_tokens = str(read1).split('\t')
            r2_tokens = str(read2).split('\t')
            r1_tags = ast.literal_eval(r1_tokens[11])
            r2_tags = ast.literal_eval(r2_tokens[11])
            r1_md = [x[1] for x in r1_tags if x[0] == 'MD'][0]
            r2_md = [x[1] for x in r2_tags if x[0] == 'MD'][0]
            r1_cigar = r1_tokens[5]
            r2_cigar = r2_tokens[5]
            r1_pos = int(r1_tokens[3])+1
            r2_pos = int(r2_tokens[3])+1
            r1_seq = r1_tokens[9]
            r2_seq = r2_tokens[9]
            # Extract data from each read
            r1_insertions, r1_deletions, r1_junctions, r1_exons, r1_mismatches = indels_junctions_exons_mismatches(r1_cigar, r1_md, r1_pos, r1_seq)
            r2_insertions, r2_deletions, r2_junctions, r2_exons, r2_mismatches = indels_junctions_exons_mismatches(r2_cigar, r2_md, r2_pos, r2_seq)
            # Process insertions
            for insertion in set(r1_insertions + r2_insertions):
                variant = (contig, insertion[0], '', insertion[1], 'I')
                if variant in all_mutations:
                    expressed_variants[variant] += 1
            # Process deletions
            for deletion in set(r1_deletions + r2_deletions):
                variant = (contig, insertion[0], deletion[1], len(deletion[1]), 'D')
                if variant in all_mutations:
                    expressed_variants[variant] += 1
            # Process mismatches
            for mismatch in set(r1_mismatches + r2_mismatches):
                ref_base = reference_index.get_stretch(contig, mismatch[0] - 1, len(mismatch[1]))
                variant = (contig, mismatch[0], ref_base, mismatch[1], 'V')
                if variant in all_mutations:
                    expressed_variants[variant] += 1
            # Find genomic intervals covered, accounting for junctions
            intervals = []
            for e in r1_exons + r2_exons:
                intervals.extend([e[0], e[1]+1])
            # Find and store variants that overlap intervals
            var_set = set()
            for i in range(0, len(intervals), 2):
                overlapping_variants = var_intervals[search_contig].overlap(intervals[i], intervals[i+1])
                for var in [x.data for x in overlapping_variants]:
                    var_set.add(var)
            # Add read counts to all overlapping variants
            for var in var_set:
                covered_variants[var] += 1
    # Remove temporary files
    if remove_files:
        for file_to_remove in files_to_remove:
            os.remove(file_to_remove)
    return expressed_variants, covered_variants



