#!/usr/bin/env python
# coding=utf-8
"""
transcript.py

Part of neoepiscope
Functions for manipulating transcripts and enumerating neoepitopes.

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
from . import bowtie_index
import collections
import copy
import bisect
import string
import re
import os
import io
import pickle
from intervaltree import IntervalTree
from operator import itemgetter
from numpy import median
import sys
import warnings
import contextlib
import networkx as nx

revcomp_translation_table = str.maketrans("ATCGI", "TAGCI")


@contextlib.contextmanager
def xopen(gzipped, *args):
    """Passes args on to the appropriate opener, gzip or regular.
    In compressed mode, functionality almost mimics gzip.open,
    but uses gzip at command line.
    As of PyPy 2.5, gzip.py appears to leak memory when writing to
    a file object created with gzip.open().
    gzipped: True iff gzip.open() should be used to open rather than
        open(); False iff open() should be used; None if input should be
        read and guessed; '-' if writing to stdout
    *args: unnamed arguments to pass
    Yield value: file object
    """
    import sys

    if gzipped == "-":
        fh = sys.stdout
    else:
        if not args:
            raise IOError("Must provide filename")
        import gzip

        if gzipped is None:
            with open(args[0], "rb") as binary_input_stream:
                # Check for magic number
                b2 = binary_input_stream.read(2)
                if b2 == b"\x1f\x8b":
                    gzipped = True
                else:
                    gzipped = False
        if gzipped:
            try:
                mode = args[1]
            except IndexError:
                mode = "r"
            if "r" in mode:
                # Be forgiving of gzips that end unexpectedly
                # old_read_eof = gzip.GzipFile._read_eof
                # gzip.GzipFile._read_eof = lambda *args, **kwargs: None
                fh = io.TextIOWrapper(gzip.open(*args, "r"))
            elif "w" in mode or "a" in mode:
                try:
                    compresslevel = int(args[2])
                except IndexError:
                    compresslevel = 9
                if "w" in mode:
                    output_stream = open(args[0], "wb")
                else:
                    output_stream = open(args[0], "ab")
                gzip_process = subprocess.Popen(
                    ["gzip", "-%d" % compresslevel],
                    bufsize=-1,
                    stdin=subprocess.PIPE,
                    stdout=output_stream,
                )
                fh = gzip_process.stdin
            else:
                raise IOError("Mode " + mode + " not supported")
        else:
            fh = open(*args)
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()
        if "gzip_process" in locals():
            gzip_process.wait()
        if "output_stream" in locals():
            output_stream.close()
        if "old_read_eof" in locals():
            gzip.GzipFile._read_eof = old_read_eof


# From https://stackoverflow.com/questions/21222506/multiple-assignments-into-a-python-dictionary
def multiassign(d, keys, values):
    d.update(zip(keys, values))


def custom_bisect_left(a, x, lo=0, hi=None, getter=0):
    """Same as bisect.bisect_left, but compares only index "getter"
    See bisect_left source for more info.
    """

    if lo < 0:
        raise ValueError("lo must be non-negative")
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid][getter] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def kmerize_peptide(peptide, min_size=8, max_size=11, editing_positions=[], ambiguous_positions=[]):
    """Obtains subsequences of a peptide.
    normal_peptide: normal peptide seq
    min_size: minimum subsequence size
    max_size: maximum subsequence size
    editing_positions: positions where edits are made
    ambiguous_positions: positions where edits are ambiguous 
    Return value: a list of tuples of peptide substrings and
        Booleans indicating whether they contain RNA-editing and ambiguous
        RNA-editing
    """
    peptide_size = len(peptide)
    return [
        item
        for sublist in [
            [(peptide[i : i + size], True if [x for x in editing_positions if x in range(i, i+size)] else False, True if [x for x in ambiguous_positions if x in range(i, i+size)] else False) for i in range(peptide_size - size + 1)]
            for size in range(min_size, max_size + 1)
        ]
        for item in sublist
    ]

# X below denotes a stop codon
_codon_table = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "X",
    "TAG": "X",
    "TGT": "C",
    "TGC": "C",
    "TGA": "X",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


_mitochondrial_codon_table = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "X",
    "TAG": "X",
    "TGT": "C",
    "TGC": "C",
    "TGA": "W",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "M",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "X",
    "AGG": "X",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

_ambiguous_codons = ["IAA", "ICC", "IAC", "ICA", "IAI", "ICI", "IIA"]

def seq_to_peptide(
    seq,
    reverse_strand=False,
    require_ATG=False,
    mitochondrial=False,
    allow_partial_codons=False,
    return_positions=False
):
    """
    Translates nucleotide sequence into peptide sequence.
    All codons including and after stop codon are recorded as X's.
    seq: nucleotide sequence
    reverse_strand: True iff strand is -
    require_ATG: True iff search for start codon (ATG)
    mitochondrial: True iff use mitochondrial codon table
    allow_partial_codons: True iff attempt to translate partial
        codons at end of sequence
    return_positions: If True, return peptide string, editing_positions, and ambiguous_positions.

    Return value: peptide string, list of peptide warnings
    """
    if reverse_strand:
        seq = seq[::-1].translate(_complement_table)
    if require_ATG:
        start = seq.find("ATG")
        if start >= 0:
            seq = seq[start:]
        else:
            return ""
    peptide = []
    editing_positions = []
    ambiguous_positions = []
    peptide_warnings = []
    seq_size = len(seq)
    for i in range(0, seq_size - seq_size % 3, 3):
        chunk = seq[i : i + 3]
        protein_pos = round(i/3, 0)
        if 'N' not in chunk and "I" not in chunk:
            if not mitochondrial:
                codon = _codon_table[chunk]
            else:
                codon = _mitochondrial_codon_table[chunk]
        elif "I" in chunk and "N" not in chunk:
            editing_positions.append(protein_pos)
            if chunk in _ambiguous_codons:
                ambiguous_positions.append(protein_pos)
            if not mitochondrial:
                codon = _codon_table[chunk.replace("I", "G")]
            else:
                codon = _mitochondrial_codon_table[chunk.replace("I", "G")]
        elif chunk.count('N') == 1 and seq[i+2] == 'N':
             # Only 1 N in the wobble position
            if "I" not in chunk:
                # No editing 
                if not mitochondrial:
                    codon_options = set(
                        [_codon_table[''.join([seq[i : i + 2], x])]
                                                for x in 'ACGT']
                        )
                else:
                    codon_options = set(
                        [_mitochondrial_codon_table[''.join([seq[i : i + 2], x])]
                                                for x in 'ACGT']
                        )
            else:
                # Editing
                chunk_options = [''.join([seq[i:i+2], x])
                    for x in 'ACGT']
                if not mitochondrial:
                    codon_options = set([_codon_table[x.replace('I', 'G')]
                        for x in chunk_options])
                else:
                    codon_options = set([_mitochondrial_codon_table[x.replace('I', 'G')]
                        for x in chunk_options])
                editing_positions.append(protein_pos)
                if [x for x in chunk_options if x in _ambiguous_codons]:
                    ambiguous_positions.append(protein_pos)
            # Check whether amino acid can be assigned
            if len(codon_options) == 1:
                codon = list(codon_options)[0]
            else:
                codon = "?"
        else:
            # More than 1 N or N not in wobble position
            codon = "?"
        peptide.append(codon)
    if seq_size % 3:
        peptide_warnings.append("incomplete_CDS")
        if allow_partial_codons:
            # 1-2 nucleotides remaining
            if seq_size % 3 == 2:
                # 2 nucleotides remaining - check if amino acid can be determined
                if "I" not in seq[-2:]:
                    # No editing
                    if not mitochondrial:
                        codon_options = set(
                            [
                                _codon_table["".join([seq[-2:], x])]
                                for x in 'ACGT'
                            ]
                        )
                    else:
                        codon_options = set(
                            [
                                _mitochondrial_codon_table["".join([seq[-2:], x])]
                                for x in 'ACGT'
                            ]
                        )
                else:
                    # Editing
                    chunk_options = [''.join([seq[-2:], x])
                            for x in 'ACGT']
                    if not mitochondrial:
                        codon_options = set([_codon_table[x.replace('I', 'G')]
                            for x in chunk_options])
                    else:
                        codon_options = set([_mitochondrial_codon_table[x.replace('I', 'G')]
                            for x in chunk_options])
                    protein_pos = round((i+3)/3, 0)
                    editing_positions.append(protein_pos)
                    if [x for x in chunk_options if x in _ambiguous_codons]:
                        ambiguous_positions.append(protein_pos)
                if len(codon_options) == 1:
                    codon = list(codon_options)[0]
                else:
                    codon = "?"
            else:
                # Only 1 amino acid left - can't determine amino acid
                codon = "?"
            peptide.append(codon)
    if mitochondrial and peptide[0] != "M":
        peptide[0] = "M"
    if return_positions:
        return "".join(peptide), editing_positions, ambiguous_positions, peptide_warnings
    else:
        return "".join(peptide), peptide_warnings


class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype. """

    # Should we handle somatic deletions that overlap germline mutations?
    # I.E., should we break up a somatic deletion into two separate mutations
    # that surround the germline mutation? Or do we only call the somatic?

    def __init__(self, bowtie_reference_index, cds, transcript_id, seleno, rna_edit_dict=None):
        """ Initializes Transcript object.
            This class assumes edits added to a transcript are properly
            phased, consistent, and nonredundant. Most conspicuously, there
            shouldn't be SNVs or insertions among deleted bases.
            bowtie_reference_index: BowtieIndexReference object for retrieving
                reference genome sequence
            cds: list of all CDS lines for exactly one transcript from GTF;
                a line can be a list pre-split by '\t' or not yet split
            transcript_id: transcript ID
            seleno: whether translated transcript contains a selenocysteine (boolean)
            rna_edit_dict: dictionary linking transcript IDs to RNA editing sites 
        """
        assert len(cds) > 0
        if rna_edit_dict is None:
            self.rna_edit_sites = []
        else:
            try:
                self.rna_edit_sites = rna_edit_dict[transcript_id]
            except KeyError:
                self.rna_edit_sites = []
        self.bowtie_reference_index = bowtie_reference_index
        self.transcript_id = transcript_id
        self.seleno = seleno
        self.intervals = []
        self.all_transcript_warnings = []
        # Internal representation is 0-based
        self._start_codon, self._stop_codon = None, None
        # Public representation is 1-based
        self.start_codon, self.stop_codon = None, None
        last_chrom, last_strand = None, None
        for line in cds:
            if type(line) is str:
                line = line.strip().split("\t")
            try:
                assert last_chrom == line[0]
            except AssertionError:
                if last_chrom is None:
                    pass
                else:
                    raise
            try:
                assert last_strand == line[6]
            except AssertionError:
                if last_strand is None:
                    pass
                else:
                    raise
            # Use exclusive start, inclusive end 0-based coordinates internally
            if line[2] == "exon":
                self.intervals.extend([int(line[3]) - 2, int(line[4]) - 1])
            elif line[2].startswith("start_codon"):
                self.start_codon = int(line[3])
                self._start_codon = self.start_codon - 1
            elif line[2].startswith("stop_codon"):
                self.stop_codon = int(line[3])
                self._stop_codon = self.stop_codon - 1
            else:
                raise NotImplementedError("GTF sequence type not currently supported")
            last_chrom, last_strand = line[0], line[6]
        # Store edits to coding sequence only
        self.edits = collections.defaultdict(list)
        self.deletion_intervals = []
        self.chrom = last_chrom
        # Determine whether transcript on mitochondrial DNA
        if self.chrom not in ["chrM", "M", "chrMT", "MT"]:
            self.mitochondrial = False
        else:
            self.mitochondrial = True
        # Flag to indicate if there is a deletion that spans intron-exon boundary
        self.boundary_spanning_deletion = False
        self.rev_strand = True if last_strand == "-" else False
        """Assume intervals are nonoverlapping! Uncomment following lines to
        check (slower)."""
        # for i in range(1, len(self.intervals)):
        #    if self.intervals[i-1] <= self.intervals[i]:
        #        raise RuntimeError(
        #                ('CDS intervals list '
        #                 '"{}" has overlapping intervals.').format(
        #                            self.intervals
        #                        )
        #            )
        # For retrieving save point
        self.last_edits = collections.defaultdict(list)
        self.last_deletion_intervals = []
        # Need to sort to bisect_left properly when editing!
        self.intervals.sort()
        # +1 is + strand, -1 is - strand
        strand = 1 - self.rev_strand * 2
        # Get start codon index
        if self._start_codon:
            self.start_codon_index = bisect.bisect_left(
                self.intervals, self._start_codon
            )
        else:
            self.start_codon_index = None
        if self.start_codon is not None:
            # Find start codon sequence
            if (
                bisect.bisect_left(self.intervals, self.start_codon + 1)
                == self.start_codon_index
            ):
                if not self.start_codon_index % 2:
                    # Move faux start codon from intron to exon
                    if self.rev_strand:
                        # Move it completely into previous exon
                        self.start_codon = (
                            self.intervals[self.start_codon_index - 1] - 1
                        )
                    else:
                        # Move it completely into next exon
                        self.start_codon = self.intervals[self.start_codon_index] + 2
                    self._start_codon = self.start_codon - 1
                    self.start_codon_index = bisect.bisect_left(
                        self.intervals, self._start_codon
                    )
                # Entire start codon contained within same exon
                self.start_codon_seq = self.bowtie_reference_index.get_stretch(
                    self.chrom, self._start_codon, 3
                )
                self.start_coordinates = [
                    x for x in range(self.start_codon, self.start_codon + 3)
                ]
            else:
                # Start codon is split across exons
                if (
                    bisect.bisect_left(self.intervals, self.start_codon)
                    == self.start_codon_index
                ):
                    # Only need to grab 1 nucleotide from next exon
                    self.start_codon_seq = self.bowtie_reference_index.get_stretch(
                        self.chrom, self._start_codon, 2
                    )
                    self.start_codon_seq += self.bowtie_reference_index.get_stretch(
                        self.chrom, self.intervals[self.start_codon_index + 1] + 1, 1
                    )
                    self.start_coordinates = [
                        self.start_codon,
                        self.start_codon + 1,
                        self.intervals[self.start_codon_index + 1] + 2,
                    ]
                else:
                    # Need to grab two nucleotides from following exon(s)
                    self.start_codon_seq = self.bowtie_reference_index.get_stretch(
                        self.chrom, self._start_codon, 1
                    )
                    self.start_coordinates = [self.start_codon]
                    if (
                        self.intervals[self.start_codon_index + 2]
                        - self.intervals[self.start_codon_index + 1]
                        >= 2
                    ):
                        # The remainder of the start codon is within the next exon
                        self.start_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom,
                            self.intervals[self.start_codon_index + 1] + 1,
                            2,
                        )
                        self.start_coordinates.extend(
                            [
                                self.intervals[self.start_codon_index + 1] + 2,
                                self.intervals[self.start_codon_index + 1] + 3,
                            ]
                        )
                    else:
                        # The start codon is split across three exons
                        self.start_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom,
                            self.intervals[self.start_codon_index + 1] + 1,
                            1,
                        )
                        self.start_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom,
                            self.intervals[self.start_codon_index + 3] + 1,
                            1,
                        )
                        self.start_coordinates.extend(
                            [
                                self.intervals[self.start_codon_index + 1] + 2,
                                self.intervals[self.start_codon_index + 3] + 2,
                            ]
                        )
            if self.rev_strand:
                # Reverse strand transcript - reverse complement the sequence
                self.start_codon_seq = self.start_codon_seq[::-1].translate(
                    revcomp_translation_table
                )
        else:
            # There is no start codon sequence
            self.start_codon_seq = None
            self.start_coordinates = []
        # Get stop codon index
        if self.stop_codon:
            self.stop_codon_index = bisect.bisect_left(self.intervals, self._stop_codon)
        else:
            self.stop_codon_index = None
        # self.stop_past_tx = False
        if self.stop_codon is not None:
            # Find stop codon sequence
            if (
                bisect.bisect_left(self.intervals, self.stop_codon + 1)
                == self.stop_codon_index
            ):
                # Entire stop codon contained within same exon
                self.stop_codon_seq = self.bowtie_reference_index.get_stretch(
                    self.chrom, self._stop_codon, 3
                )
                self.stop_coordinates = [
                    x for x in range(self.stop_codon, self.stop_codon + 3)
                ]
            else:
                # Stop codon is split across exons
                if (
                    bisect.bisect_left(self.intervals, self.stop_codon)
                    == self.stop_codon_index
                ):
                    # Only need to grab 1 nucleotide from following exon
                    self.stop_codon_seq = self.bowtie_reference_index.get_stretch(
                        self.chrom, self._stop_codon, 2
                    )
                    self.stop_codon_seq += self.bowtie_reference_index.get_stretch(
                        self.chrom, self.intervals[self.stop_codon_index + 1] + 1, 1
                    )
                    self.stop_coordinates = [
                        self.stop_codon,
                        self.stop_codon + 1,
                        self.intervals[self.stop_codon_index + 1] + 2,
                    ]
                else:
                    # Need to grab two nucleotides from following exon(s)
                    self.stop_codon_seq = self.bowtie_reference_index.get_stretch(
                        self.chrom, self._stop_codon, 1
                    )
                    self.stop_coordinates = [self.stop_codon]
                    if (
                        self.intervals[self.stop_codon_index + 2]
                        - self.intervals[self.stop_codon_index + 1]
                        >= 2
                    ):
                        # The remainder of the stop codon is within the next exon
                        self.stop_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom, self.intervals[self.stop_codon_index + 1] + 1, 2
                        )
                        self.stop_coordinates.extend(
                            [
                                self.intervals[self.stop_codon_index + 1] + 2,
                                self.intervals[self.stop_codon_index + 1] + 3,
                            ]
                        )
                    else:
                        # The stop codon is split across three exons
                        self.stop_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom, self.intervals[self.stop_codon_index + 1] + 1, 1
                        )
                        self.stop_codon_seq += self.bowtie_reference_index.get_stretch(
                            self.chrom, self.intervals[self.stop_codon_index + 3] + 1, 1
                        )
                        self.stop_coordinates.extend(
                            [
                                self.intervals[self.stop_codon_index + 1] + 2,
                                self.intervals[self.stop_codon_index + 3] + 2,
                            ]
                        )
            if self.rev_strand:
                # Reverse strand transcript - reverse complement the sequence
                self.stop_codon_seq = self.stop_codon_seq[::-1].translate(
                    revcomp_translation_table
                )
        else:
            # There is no stop codon sequence
            self.stop_codon_seq = None
            self.stop_coordinates = []

    def reset(self, reference=False):
        """Resets to last save point or reference (i.e., removes all edits).
        reference: if False, tries to reset to last save point, and if that
            doesn't exist, resets to reference. If True, resets to
            reference.
        No return value.
        """
        if reference:
            self.edits = collections.defaultdict(list)
            self.deletion_intervals = []
            self.boundary_spanning_deletion = False
        else:
            self.edits = copy.deepcopy(self.last_edits)
            self.deletion_intervals = copy.deepcopy(self.last_deletion_intervals)
            if [x for x in self.deletion_intervals if x[3][6]] == []:
                self.boundary_spanning_deletion = False

    def edit(self, seq, pos, mutation_type="V", mutation_class="S", vaf=None):
        """ Adds an edit to the transcript.
            seq: sequence to add or delete from reference; for deletions, all
                that matters is this sequence has the same length as the
                sequence to delete. Also for deletions, seq can be an integer
                specifying how many bases to delete.
            pos: 1-based coordinate. For insertions, this is the coordinate
                directly before the inserted sequence. For deletions, this
                is the coordinate of the first base of the transcript to be
                deleted. Coordinates are always w.r.t. genome.
            mutation_type: V for SNV, I for insertion, D for deletion,
                           E for RNA edit.
            mutation_class: S for somatic, G for germline, P for post-transcriptional
            vaf: variant allele frequency (None if not available)
            No return value.
        """
        ## Need to add check for only 1 mutation of each class per position
        if mutation_type == "D":
            try:
                deletion_size = int(seq)
            except ValueError:
                deletion_size = len(seq)
                ref_deletion = self.bowtie_reference_index.get_stretch(
                    self.chrom, pos - 1, pos + deletion_size + 1 - pos - 1
                )
                if (bisect.bisect_left(self.intervals, pos - 2) % 2) and (
                    bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                ):
                    boundary_span = False
                elif (not (bisect.bisect_left(self.intervals, pos - 2) % 2)) and (
                    not (
                        bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                    )
                ):
                    boundary_span = False
                elif (pos - 2) in self.intervals and (
                    bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                ):
                    boundary_span = False
                elif (pos + deletion_size - 2) in self.intervals and (
                    bisect.bisect_left(self.intervals, pos - 2) % 2
                ):
                    boundary_span = False
                else:
                    boundary_span = True
                    self.boundary_spanning_deletion = True
                if seq == ref_deletion:
                    del_interval = (
                        pos - 2,
                        pos + deletion_size - 2,
                        mutation_class,
                        (self.chrom, pos, seq, "", mutation_type, vaf, boundary_span),
                    )
                else:
                    raise RuntimeError(
                        "".join(
                            [
                                "Deletion of ",
                                seq,
                                " at position ",
                                str(pos),
                                " on contig ",
                                self.chrom,
                                " is incompatible with reference",
                            ]
                        )
                    )
            else:
                if (bisect.bisect_left(self.intervals, pos - 2) % 2) and (
                    bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                ):
                    boundary_span = False
                elif (not (bisect.bisect_left(self.intervals, pos - 2) % 2)) and (
                    not (
                        bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                    )
                ):
                    boundary_span = False
                elif (pos - 2) in self.intervals and (
                    bisect.bisect_left(self.intervals, pos + deletion_size - 2) % 2
                ):
                    boundary_span = False
                elif (pos + deletion_size - 2) in self.intervals and (
                    bisect.bisect_left(self.intervals, pos - 2) % 2
                ):
                    boundary_span = False
                else:
                    boundary_span = True
                    self.boundary_spanning_deletion = True
                del_interval = (
                    pos - 2,
                    pos + deletion_size - 2,
                    mutation_class,
                    (
                        self.chrom,
                        pos,
                        self.bowtie_reference_index.get_stretch(
                            self.chrom, pos - 1, pos + deletion_size + 1 - pos - 1
                        ),
                        "",
                        mutation_type,
                        vaf,
                        boundary_span,
                    ),
                )
            for interval in self.deletion_intervals:
                if del_interval[2] == interval[2]:
                    if (
                        del_interval[0] >= interval[0] and del_interval[0] < interval[1]
                    ) or (
                        del_interval[1] > interval[0] and del_interval[1] <= interval[1]
                    ):
                        class_dict = {"S": "somatic", "G": "germline"}
                        raise NotImplementedError(
                            "".join(
                                [
                                    "2 deletions of same class cannot",
                                    " overlap - was deletion of ",
                                    del_interval[3][2],
                                    " at ",
                                    str(pos),
                                    " or deletion of ",
                                    interval[3][2],
                                    " at ",
                                    str(interval[0] + 2),
                                    " not a ",
                                    class_dict[mutation_class],
                                    " mutation?",
                                ]
                            )
                        )
            self.deletion_intervals.append(del_interval)
        elif mutation_type == "I":
            other_insertions = [edit for edit in self.edits[pos - 1] if edit[1] == "I"]
            if len(other_insertions) == 0:
                self.edits[pos - 1].append(
                    (
                        seq,
                        mutation_type,
                        mutation_class,
                        (self.chrom, pos, "", seq, mutation_type, vaf),
                    )
                )
            else:
                raise NotImplementedError(
                    "".join(
                        [
                            "2 insertions cannot be added at same position",
                            " - is insertion of ",
                            seq,
                            " at ",
                            str(pos),
                            " valid?",
                        ]
                    )
                )
        elif mutation_type == "V":
            reference_seq = self.bowtie_reference_index.get_stretch(
                self.chrom, pos - 1, len(seq)
            )
            other_snvs = [edit for edit in self.edits[pos - 1] if edit[1] == "V"]
            if mutation_class not in [snv[2] for snv in other_snvs]:
                self.edits[pos - 1].append(
                    (
                        seq,
                        mutation_type,
                        mutation_class,
                        (self.chrom, pos, reference_seq, seq, mutation_type, vaf),
                    )
                )
            else:
                class_dict = {"S": "somatic", "G": "germline"}
                raise NotImplementedError(
                    "".join(
                        [
                            "2 SNVs of same class cannot ",
                            "be added at same position",
                            " - was mutation of ",
                            reference_seq,
                            " to ",
                            seq,
                            " at ",
                            str(pos),
                            " not a ",
                            class_dict[mutation_class],
                            " mutation?",
                        ]
                    )
                )
        elif mutation_type == "E":
            reference_seq = self.bowtie_reference_index.get_stretch(
                self.chrom, pos - 1, len(seq)
            )
            existing_RNA_edit = [edit for edit in self.edits[pos - 1] if edit[1] == 'E']
            if not existing_RNA_edit:
                if not self.rev_strand and self.start_coordinates[0] == pos:
                    self.all_transcript_warnings.append("rna_editing_may_disrupt_start_codon")
                    warnings.warn("Start codon in %s may be disrupted by RNA editing" % self.transcript_id)
                elif self.rev_strand and self.start_coordinates[2] == pos:
                    self.all_transcript_warnings.append("rna_editing_may_disrupt_start_codon")
                    warnings.warn("Start codon in %s may be disrupted by RNA editing" % self.transcript_id)
                self.edits[pos - 1].append(
                    (
                        seq,
                        mutation_type,
                        mutation_class,
                        (self.chrom, pos, reference_seq, seq, mutation_type, vaf),
                    )
                )
            else:
                warnings.warn("Cannot add RNA edit to a site in %s that already has an RNA edit" 
                    % self.transcript_id)
        else:
            raise NotImplementedError("Mutation type not yet implemented")

    def expressed_edits(
        self, start=None, end=None, genome=True, include_somatic=1, include_germline=2, 
        include_rna_edits=0):
        """ Gets expressed set of edits and transcript intervals.
            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified
            include_somatic: 0 = do not include somatic mutations,
                1 = exclude somatic mutations from reference comparison,
                2 = include somatic mutations in both annotated sequence and
                reference comparison
            include_germline: 0 = do not include germline mutations,
                1 = exclude germline mutations from reference comparison,
                2 = include germline mutations in both annotated sequence and
                reference comparison
            include_rna_edits: 0 = do not include A to I RNA editing,
                1 = exclude RNA edits from reference comparison,
                2 = include RNA edits in both annotated sequence and reference comparison
            Return value: tuple (defaultdict
                                 mapping edits to lists of
                                 (seq, mutation_type, mutation_class)
                                 tuples, interval list; this is a list of
                                 tuples (bound, {'R', 'G', or 'S'}), which
                                 says whether the bound is due to CDS bound
                                 ("R"), a germline deletion ("G"), or a
                                 somatic deletion ("S"))
        """
        if not genome:
            raise NotImplementedError(
                "Retrieving sequence with transcript coordinates not "
                "yet fully supported."
            )
        if start is None:
            start = self.intervals[0] + 1
        else:
            start -= 1
        if end is None:
            end = self.intervals[-1]
        else:
            end -= 1
        assert end >= start
        self.save()
        if include_rna_edits:
            for item in self.rna_edit_sites:
                self.edit("I", item[1], mutation_type = "E", mutation_class="P")
        # Change start and end intervals of CDS intervals
        start_index = bisect.bisect_left(self.intervals, start)
        if not (start_index % 2):
            # start should be beginning of a CDS
            start_index += 1
            try:
                start = self.intervals[start_index - 1] + 1
            except IndexError:
                # Start is outside bounds of transcript
                return ""
        end_index = bisect.bisect_left(self.intervals, end)
        if not (end_index % 2):
            # end should be end of CDS
            end = self.intervals[end_index - 1]
            end_index -= 1
        intervals = [start - 1] + self.intervals[start_index:end_index] + [end]
        assert len(intervals) % 2 == 0
        # Include only relevant deletion intervals
        relevant_deletion_intervals, edits = [], collections.defaultdict(list)
        for i in range(0, len(self.deletion_intervals)):
            if (
                self.deletion_intervals[i][0] <= self.intervals[0]
                and self.deletion_intervals[i][1] >= self.intervals[-1]
            ):
                self.reset()
                return (
                    {},
                    [
                        (self.intervals[0], "R", ()),
                        (
                            self.intervals[0],
                            self.deletion_intervals[i][2],
                            self.deletion_intervals[i][3],
                        ),
                        (self.intervals[-1], "R", ()),
                        (self.intervals[-1], "R", ()),
                    ],
                )
        sorted_deletion_intervals = [
            interval
            for interval in self.deletion_intervals
            if (
                interval[2] == "S"
                and include_somatic
                or interval[2] == "G"
                and include_germline
            )
        ]
        if sorted_deletion_intervals:
            sorted_deletion_intervals.sort(key=itemgetter(0, 1))
            deletion_intervals = [
                (
                    sorted_deletion_intervals[0][0],
                    sorted_deletion_intervals[0][2],
                    sorted_deletion_intervals[0][3],
                ),
                (
                    sorted_deletion_intervals[0][1],
                    sorted_deletion_intervals[0][2],
                    sorted_deletion_intervals[0][3],
                ),
            ]
            for i in range(1, len(sorted_deletion_intervals)):
                if sorted_deletion_intervals[i][0] <= deletion_intervals[-1][0]:
                    deletion_intervals[-2] = min(
                        deletion_intervals[-2],
                        (
                            sorted_deletion_intervals[i][0],
                            sorted_deletion_intervals[i][2],
                            sorted_deletion_intervals[i][3],
                        ),
                        key=itemgetter(0),
                    )
                    deletion_intervals[-1] = max(
                        deletion_intervals[-1],
                        (
                            sorted_deletion_intervals[i][1],
                            sorted_deletion_intervals[i][2],
                            sorted_deletion_intervals[i][3],
                        ),
                        key=itemgetter(0),
                    )
                else:
                    deletion_intervals.extend(
                        [
                            (
                                sorted_deletion_intervals[i][0],
                                sorted_deletion_intervals[i][2],
                                sorted_deletion_intervals[i][3],
                            ),
                            (
                                sorted_deletion_intervals[i][1],
                                sorted_deletion_intervals[i][2],
                                sorted_deletion_intervals[i][3],
                            ),
                        ]
                    )
            for i in range(0, len(deletion_intervals), 2):
                start_index = bisect.bisect_left(intervals, deletion_intervals[i][0])
                end_index = bisect.bisect_left(intervals, deletion_intervals[i + 1][0])
                if start_index == end_index:
                    if start_index % 2:
                        # Entirely in a single interval
                        relevant_deletion_intervals.extend(
                            deletion_intervals[i : i + 2]
                        )
                    # else deletion is entirely outside CDS within start/end
                elif (
                    start_index == (end_index - 1)
                    and deletion_intervals[i][0] == intervals[start_index]
                ):
                    relevant_deletion_intervals.extend(deletion_intervals[i : i + 2])
                elif start_index == 0:
                    relevant_deletion_intervals.append((intervals[0], "R", tuple()))
                    if end_index % 2:
                        end_pos = deletion_intervals[i + 1]
                        relevant_deletion_intervals.extend(
                            [
                                (intervals[i], "R", tuple())
                                for i in range(start_index + 1, end_index)
                            ]
                        )
                    else:
                        end_pos = (intervals[end_index - 1], "R", tuple())
                        relevant_deletion_intervals.extend(
                            [
                                (intervals[i], "R", tuple())
                                for i in range(start_index, end_index)
                            ]
                        )
                    relevant_deletion_intervals.append(end_pos)
                else:
                    assert end_index > start_index
                    if (
                        not start_index % 2
                        and end_index % 2
                        and start_index == (end_index - 1)
                    ):
                        relevant_deletion_intervals.extend(
                            [
                                (intervals[start_index], "R", tuple()),
                                (
                                    deletion_intervals[i + 1][0] - 1,
                                    deletion_intervals[i + 1][1],
                                    deletion_intervals[i + 1][2],
                                ),
                            ]
                        )
                    elif not start_index % 2 and not end_index % 2:
                        relevant_deletion_intervals.extend(
                            [
                                (intervals[i], "R", tuple())
                                for i in range(start_index, end_index)
                            ]
                        )
                    else:
                        if (
                            start_index
                            % 2
                            # or deletion_intervals[i][0] == intervals[start_index]
                        ):
                            pos = deletion_intervals[i]
                        else:
                            pos = (intervals[start_index], "R", tuple())
                            start_index += 1
                        # deletion_intervals[i] becomes a new end
                        relevant_deletion_intervals.extend(
                            [pos, (intervals[start_index], "R", tuple())]
                        )
                        if end_index % 2:
                            end_pos = deletion_intervals[i + 1]
                            relevant_deletion_intervals.extend(
                                [
                                    (intervals[i], "R", tuple())
                                    for i in range(start_index + 1, end_index)
                                ]
                            )
                        else:
                            end_pos = (intervals[end_index - 1], "R", tuple())
                            relevant_deletion_intervals.extend(
                                [
                                    (intervals[i], "R", tuple())
                                    for i in range(start_index, end_index)
                                ]
                            )
                        relevant_deletion_intervals.append(end_pos)
        intervals = sorted(
            [(interval, "R", tuple()) for interval in intervals]
            + relevant_deletion_intervals,
            key=lambda k: (k[0], k[1] != "R"),
        )
        # Remove empty intervals
        intervals = [
            intervals[i]
            for i in range(len(intervals))
            if (
                i % 2
                and intervals[i][0] != intervals[i - 1][0]
                or i % 2 == 0
                and intervals[i + 1][0] != intervals[i][0]
            )
        ]
        # Only associate one end of a deletion interval with deletion
        #   to prevent including it multiple times
        adjusted_intervals = [intervals[0]]
        deletion_data = []
        if intervals[0][1] != "R":
            deletion_data.append(intervals[0][2])
        for i in range(1, len(intervals)):
            if intervals[i][1] == "R":
                adjusted_intervals.append(intervals[i])
            else:
                if intervals[i][2] not in deletion_data:
                    adjusted_intervals.append(intervals[i])
                    deletion_data.append(intervals[i][2])
                else:
                    adjusted_intervals.append((intervals[i][0], "R", tuple(), None))
        edits = collections.defaultdict(list)
        for pos in self.edits:
            # Add edit if and only if it's in one of the CDSes
            start_index = custom_bisect_left(intervals, pos)
            for edit in self.edits[pos]:
                if (
                    include_somatic
                    and edit[2] == "S" # somatic
                    or include_germline
                    and edit[2] == "G" # germline
                ):
                    if edit[1] == "V": # substitution
                        ## need to determine why 'and' clause part of if statement
                        if start_index % 2 and edit[3][1] != edit[0]:
                            # Add edit if and only if it lies within bounds
                            edits[pos].append(edit)
                    elif edit[1] == "I": # insertion
                        try:
                            if start_index % 2 or pos == intervals[start_index][0]:
                                # An insertion is valid before or after a block
                                edits[pos].append(edit)
                        except IndexError:
                            continue
                elif (
                    include_rna_edits
                    and edit[1] == "E" # RNA edit
                ):
                    if start_index % 2 and edit[3][1] != edit[0]:
                        edits[pos].append(edit)
        # Handle germline, somatic variants and RNA edits (if desired) at same pos
        edits_to_return = copy.deepcopy(edits)
        for pos, edits_at_pos in edits.items():
            edits_at_pos = [x for x in edits_at_pos if x[1] in "VE"]
            # Two or more overlapping variants, eg. germline + somatic SNVs, germline SNV + RNA edit
            if len(edits_at_pos) > 1:
                edits_at_pos = sorted(edits[pos], key=lambda x: (x[1], x[2]))
                new_entry = [x for x in edits[pos] if x[1] == "I"]
                if edits_at_pos[0][1] == 'V':
                    germline = edits_at_pos[0]
                    somatic = edits_at_pos[1]
                    mutation_type = "V"
                    seq = somatic[0]
                    mutation_class = "S"
                    var = list(somatic[3])
                    if include_somatic == 2:
                        var[2] = somatic[3][3]
                    elif include_germline == 2:
                        #include_somatic == 1
                        var[2] = germline[3][3]
                else:
                    # RNA edits present
                    germline = (edits_at_pos[1] if edits_at_pos[1][2] == 'G'
                                else None)
                    if germline is None:
                        somatic = edits_at_pos[1]
                    elif len(edits_at_pos) > 2:
                        somatic = edits_at_pos[2]
                    else:
                        somatic = None
                    mutation_type = "V"
                    if germline and somatic:
                        seq = somatic[0]
                        mutation_class = "S"
                        var = list(somatic[3])
                        if include_somatic == 2:
                            var[2] = somatic[3][3]
                        elif include_germline == 2:
                            # include_somatic == 1
                            var[2] = germline[3][3]
                    elif germline:
                        # include_germline == 1
                        seq = germline[0]
                        mutation_class = "G"
                        var = list(germline[3])
                        if include_germline == 2:
                            var[2] = germline[3][3]
                    else:
                        # include_somatic == 1
                        seq = somatic[0]
                        mutation_class = "S"
                        var = list(somatic[3])
                        if include_somatic == 2:
                            var[2] = somatic[3][3]
                    # current behavior: everything gets RNA-edited if include_rna_edits
                    # future behavior: include_rna_edits = [0,1,2] where behavior would change based on specification of where to include edits
                    if include_rna_edits == 1:
                        if (seq == 'T' and self.rev_strand or
                            seq == 'A' and not self.rev_strand):
                            # Favor RNA edit leaving ref as-is
                            seq = 'I'
                            mutation_type = "E"
                            var[3] = seq
                            var[4] = "E"
                    if include_rna_edits == 2:
                        if (seq == 'T' and self.rev_strand or
                            seq == 'A' and not self.rev_strand):
                            # Favor RNA edit
                            seq = 'I'
                            mutation_type = "E"
                            var[3] = seq
                            var[4] = "E"
                        if (var[2] == 'T' and self.rev_strand or
                            var[2] == 'A' and not self.rev_strand):
                            # RNA edit becomes ref
                            var[2] = 'I'
                # Update entry
                new_entry.append(
                        (seq, mutation_type, mutation_class, tuple(var))
                    )
                edits_to_return[pos] = new_entry
            # RNA-edit present
            elif len(edits_at_pos) == 1 and edits_at_pos[0][1] == "E":
                if not include_rna_edits:
                    del edits_to_return[pos]
                else:
                    ref_at_pos = edits_at_pos[0][3][2]
                    if not (ref_at_pos == 'T' and self.rev_strand or
                            ref_at_pos == 'A' and not self.rev_strand):
                        warnings.warn("Reference nucleotide is not A or T at RNA "
                                      "edit site {}:{}; ignoring.".format(
                                                self.chrom, pos
                                            )
                                        )
                        del edits_to_return[pos]
                    elif include_rna_edits == 1:
                        edits_to_return[pos] = [
                            ('I', 'E', 'P', (el[3][0], el[3][1], el[3][2], 'I', 'E'))
                                for el in edits_to_return[pos]
                        ]
                    else:
                        edits_to_return[pos] = [
                            ('I', 'E', 'P', (el[3][0], el[3][1], 'I', 'I', 'E'))
                                for el in edits_to_return[pos]
                        ]
        self.reset()
        return (edits_to_return, adjusted_intervals)

    def save(self):
        """Creates save point for edits.

        No return value.
        """
        self.last_edits = copy.deepcopy(self.edits)
        self.last_deletion_intervals = copy.deepcopy(self.deletion_intervals)

    def reading_frame(self, pos):
        """Retrieves reading frame (0, 1, or 2) at given coordinate.

        NOTE: must be updated to include chromosome to accommodate fusions
        pos: 1-based position at which reading frame is desired
        Return value: reading frame; 0 means first base of codon, 1 means
        second base, and 2 means third base. None means the coordinate is
        outside the coding sequence of a given transcript.
        """
        pos -= 1
        pos_index = bisect.bisect_left(self.intervals, pos)
        if (
            not (pos_index % 2)
            or not self.start_codon_index
            # or not self.stop_codon_index
        ):
            # We're outside exon sequence
            return None
        if self.rev_strand:
            if pos_index == self.start_codon_index:
                # Within the same interval as the start codon
                if pos > self._start_codon + 2:
                    # Outside coding sequence
                    return None
                return (self._start_codon + 2 - pos) % 3
            else:
                if pos > self._start_codon or (
                    self.stop_codon is not None and pos < self._stop_codon
                ):
                    return None
                seq_length = (
                    (self.intervals[pos_index] - pos + 1)
                    + (
                        self._start_codon
                        + 2
                        - self.intervals[self.start_codon_index - 1]
                    )
                    + sum(
                        [
                            self.intervals[i + 1] - self.intervals[i]
                            for i in range(pos_index + 1, self.start_codon_index - 1, 2)
                        ]
                    )
                )
                return (seq_length - 1) % 3
        else:
            if pos_index == self.start_codon_index:
                if pos < self._start_codon:
                    return None
                return (pos - self._start_codon) % 3
            else:
                if pos < self._start_codon or (
                    self.stop_codon is not None and pos > self._stop_codon
                ):
                    return None
                seq_length = (
                    (pos - self.intervals[pos_index - 1])
                    + (self.intervals[self.start_codon_index] - self._start_codon + 1)
                    + sum(
                        [
                            self.intervals[i + 1] - self.intervals[i]
                            for i in range(self.start_codon_index + 1, pos_index - 1, 2)
                        ]
                    )
                )
                return (seq_length - 1) % 3

    def _seq_append(
        self, seq_list, seq, mutation_class, mutation_info, position, merge=True
    ):
        """Appends mutation to seq_list, merging successive mutations.
        seq_list: list of tuples (sequence, type) where type is one
            of R, G, S, or E (for respectively reference, germline edit,
            somatic edit, or RNA edit). Empty sequence means there was a deletion.
        seq: seq to add
        mutation_class: S for somatic, G for germline, R for reference, P for post-transcriptional
        mutation_info: tuple containing (1 based mutation position from
            vcf, mutation sequence, mutation type, and VAF)
        position: 1-based genomic position of first base added
        No return value; seq_list is merely updated.
        """
        if not merge:
            if seq or mutation_class != "R":
                if isinstance(mutation_info, list):
                    seq_list.append(
                        (
                            seq,
                            mutation_class,
                            [mutation_info[i] for i in range(0, len(mutation_info))],
                            position,
                        )
                    )
                else:
                    seq_list.append((seq, mutation_class, [mutation_info], position))
            return
        try:
            condition = seq_list[-1][1] == mutation_class
        except IndexError:
            # Add first item in seq_list
            assert not seq_list
            if seq or mutation_class != "R":
                if isinstance(mutation_info, list):
                    seq_list.append(
                        (
                            seq,
                            mutation_class,
                            [mutation_info[i] for i in range(0, len(mutation_info))],
                            position,
                        )
                    )
                else:
                    seq_list.append((seq, mutation_class, [mutation_info], position))
            return
        if condition:
            if self.rev_strand:
                adjusted_position = position
            else:
                adjusted_position = seq_list[-1][3]
            if isinstance(mutation_info, list):
                adjusted_mutation_info = seq_list[-1][2] + [
                    x for x in mutation_info if x not in seq_list[-1][2]
                ]
            elif mutation_info not in seq_list[-1][2]:
                adjusted_mutation_info = seq_list[-1][2] + [mutation_info]
            else:
                adjusted_mutation_info = seq_list[-1][2]
            seq_list[-1] = (
                seq_list[-1][0] + seq,
                mutation_class,
                adjusted_mutation_info,
                adjusted_position,
            )
        elif seq or mutation_class != "R":
            if isinstance(mutation_info, list):
                seq_list.append(
                    (
                        seq,
                        mutation_class,
                        [mutation_info[i] for i in range(len(mutation_info))],
                        position,
                    )
                )
            else:
                seq_list.append((seq, mutation_class, [mutation_info], position))

    def hybridize_seq(self, prev_seq, new_seq, include_somatic=1, include_germline=2):
        """Hybridizes overlapping germline and somatic deletions for
        annotated_seq()

        prev_seq: the previous deletion
        new_seq: the new deletion which overlaps with prev_seq
        reverse_strand: boolean (whether or
                        not transcript is reverse strand)
        include_somatic: 0 = do not include somatic mutations,
            1 = exclude somatic mutations from reference comparison,
            2 = include somatic mutations in both annotated sequence and
            reference comparison
        include_germline: 0 = do not include germline mutations,
            1 = exclude germline mutations from reference comparison,
            2 = include germline mutations in both annotated sequence and
            reference comparison

        Return value: hybridized sequence - tuple of ('', 'H',
            mutation_information, position) where 'H' denotes a hybrid of
            somatic and germline deletions, mutation_information is a list
            of [ALT, REF], where ALT and REF are lists of [(chr, adj_pos,
            adj_seq [meta information])] for the alternate and reference
            sequences, respectively, and position is the 1-based genomic
            position of the beginning of the hybrid deletion
        """
        if type(prev_seq[2][0]) == list:
            # Combining a third or subsequent deletion with 2 or more previous
            alt, ref = prev_seq[2][0], prev_seq[2][1]
        else:
            # Creating new hybrid deletion
            assert sorted([prev_seq[1], new_seq[1]]) == ["G", "S"]
            prev_chrom, prev_pos = prev_seq[2][0][0], prev_seq[3]
            prev_del_seq = "".join([x[2] for x in prev_seq[2]])
            prev_mut_info = prev_seq[2]
            alt = [prev_chrom, prev_pos, prev_del_seq, "", prev_mut_info]
            if (prev_seq[1] == "G" and include_germline == 2) or (
                prev_seq[1] == "S" and include_somatic == 2
            ):

                ref = [self.chrom, prev_pos, prev_del_seq, "", prev_mut_info]
            else:
                ref = [self.chrom, new_seq[3], "", prev_del_seq, []]
        adj_alt_seq = copy.deepcopy(alt[2])
        adj_alt_allele = copy.deepcopy(alt[3])
        adj_alt_mut_info = copy.deepcopy(alt[4])
        adj_ref_seq = copy.deepcopy(ref[2])
        adj_ref_allele = copy.deepcopy(ref[3])
        adj_ref_mut_info = copy.deepcopy(ref[4])
        if self.rev_strand:
            for mut in new_seq[2]:
                adj_index = max(
                    i
                    for i in range(len(mut[2]) + 1)
                    if adj_alt_seq[::-1]
                    .translate(revcomp_translation_table)
                    .endswith(mut[2][::-1].translate(revcomp_translation_table)[:i])
                )
                added_seq = (
                    mut[2][::-1]
                    .translate(revcomp_translation_table)[adj_index:][::-1]
                    .translate(revcomp_translation_table)
                )
                adj_alt_seq = added_seq + adj_alt_seq
                adj_alt_mut_info.append(mut)
                alt[1] = new_seq[3]
                adj_ref_index = max(
                    i
                    for i in range(len(mut[2]) + 1)
                    if adj_ref_seq[::-1]
                    .translate(revcomp_translation_table)
                    .endswith(mut[2][::-1].translate(revcomp_translation_table)[:i])
                )
                added_ref_seq = (
                    mut[2][::-1]
                    .translate(revcomp_translation_table)[adj_index:][::-1]
                    .translate(revcomp_translation_table)
                )
                if (new_seq[1] == "G" and include_germline == 2) or (
                    new_seq[1] == "S" and include_somatic == 2
                ):
                    adj_ref_seq = added_ref_seq + adj_ref_seq
                    adj_ref_mut_info.append(mut)
                    ref[1] = new_seq[3]
                    if adj_ref_seq != "":
                        index = prev_seq[3] - new_seq[3]
                        adj_ref_allele = adj_ref_allele[:index]
                else:
                    adj_ref_allele = added_ref_seq + adj_ref_allele
        else:
            for mut in new_seq[2]:
                adj_alt_seq += mut[2][
                    max(
                        i
                        for i in range(len(mut[2]) + 1)
                        if adj_alt_seq.endswith(mut[2][:i])
                    ) :
                ]
                adj_alt_mut_info.append(mut)
                if (new_seq[1] == "G" and include_germline == 2) or (
                    new_seq[1] == "S" and include_somatic == 2
                ):
                    adj_ref_seq += mut[2][
                        max(
                            i
                            for i in range(len(mut[2]) + 1)
                            if adj_ref_seq.endswith(mut[2][:i])
                        ) :
                    ]
                    adj_ref_mut_info.append(mut)
                    if adj_ref_allele != "":
                        index = new_seq[3] - prev_seq[3]
                        adj_ref_allele = adj_ref_allele[:index]
                else:
                    adj_ref_allele += mut[2][
                        max(
                            i
                            for i in range(len(mut[2]) + 1)
                            if adj_ref_seq.endswith(mut[2][:i])
                        ) :
                    ]
        alt[2] = copy.deepcopy(adj_alt_seq)
        alt[3] = copy.deepcopy(adj_alt_allele)
        alt[4] = copy.deepcopy(adj_alt_mut_info)
        ref[2] = copy.deepcopy(adj_ref_seq)
        ref[3] = copy.deepcopy(adj_ref_allele)
        ref[4] = copy.deepcopy(adj_ref_mut_info)
        return ("", "H", [alt, ref], alt[1])

    def annotated_seq(
        self, start=None, end=None, genome=True, include_somatic=1,
        include_germline=2, include_rna_edits=0):
        """ Retrieves transcript sequence between start and end coordinates.
            Includes info on whether edits are somatic or germline and whether
            sequence is reference sequence.
            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified
            include_somatic: 0 = do not include somatic mutations,
                1 = exclude somatic mutations from reference comparison,
                2 = include somatic mutations in both annotated sequence and
                reference comparison
            include_germline: 0 = do not include germline mutations,
                1 = exclude germline mutations from reference comparison,
                2 = include germline mutations in both annotated sequence and
                reference comparison
            include_rna_edits: 0 = do not include A to I RNA editing,
                1 = exclude RNA edits from reference comparison,
                2 = include RNA edits in both annotated sequence and reference comparison
            Return value: list of tuples (sequence, mutation class,
                mutation information, position),
                where sequence is a segment of sequence of the (possibly)
                mutated transcript, mutation class is one of {'G', 'S', 'H', 'R', 'P'},
                where 'G' denotes germline, 'S' denotes somatic, 'H' denotes hybrid
                of somatic and germline, 'R' denotes reference sequence, and 'P' denotes post-transcriptional mod,
                mutation information is a list of tuple(s) (chromosome,
                1-based position of {first base of deletion, base before
                insertion, SNV, RNA edit}, reference sequence, variant sequence,
                {'D', 'I', 'V', 'E'}, VAF) , and position is the 1-based position
                of the first base of sequence. For hybrid, the tuple structure
                of mutation information is nested inside of a list: [[ALT], [REF]],
                where ALT and REF are structured as [chromosome, adj. position,
                adj. deletion, allele seq, [mutation information]] for the alternate and
                reference sequences, respectively.
        """
        # Use 0-based coordinates internally
        if start is None:
            start = self.intervals[0] + 2
        if end is None:
            end = self.intervals[-1] + 1
        if end < start:
            return []
        if genome:
            # Capture only sequence between start and end
            edits, intervals = self.expressed_edits(
                start,
                end,
                genome=True,
                include_somatic=include_somatic,
                include_germline=include_germline,
                include_rna_edits=include_rna_edits,
            )
            new_edits = copy.deepcopy(edits)
            """Check for insertions at beginnings of intervals, and if they're
            present, shift them to ends of previous intervals so they're
            actually added."""
            i = 0
            while i < len(intervals):
                if intervals[i][0] in edits and i:
                    assert (
                        len(edits[intervals[i][0]]) == 1
                        and edits[intervals[i][0]][0][1] == "I"
                    )
                    new_edits[intervals[i - 1][0]] = new_edits[intervals[i][0]]
                    del new_edits[intervals[i][0]]
                    """Code below would add insertion to first block,
                    but no more
                    if i do the above,
                    else:
                        intervals = [(-1, 'R', [], []),
                                     (-1, 'R', [], [])] + intervals
                        # Have to add 2 because we modified intervals above
                        i += 2
                        new_edits[-1] = new_edits[intervals[i][0]]
                        del new_edits[intervals[i][0]]"""
                i += 2
            # Grab reference sequence to pull from
            seqs = []
            for i in range(0, len(intervals), 2):
                seqs.append(
                    (
                        self.bowtie_reference_index.get_stretch(
                            self.chrom,
                            intervals[i][0] + 1,
                            intervals[i + 1][0] - intervals[i][0],
                        ),
                        (intervals[i][0] + 2, intervals[i + 1][0] + 1),
                    )
                )
            # Now build sequence in order of increasing edit position
            i = 1
            pos_group, final_seq = [], []
            for pos in sorted(new_edits.keys()) + [self.intervals[-1] + 1]:
                if pos > intervals[i][0]:
                    # Position is outside of this exon
                    last_index, last_pos = 0, intervals[i - 1][0] + 1
                    for pos_to_add in pos_group:
                        fill = pos_to_add - last_pos
                        if intervals[i - 1][1] != "R":
                            # Exonic section bounded by deletion on left - add to seq
                            if isinstance(intervals[i - 1][2], list):
                                genomic_position = min(
                                    [x[1] for x in intervals[i - 1][2]]
                                )
                            else:
                                genomic_position = intervals[i - 1][2][1]
                            self._seq_append(
                                final_seq,
                                "",
                                intervals[i - 1][1],
                                intervals[i - 1][2],
                                genomic_position,
                                merge=False,
                            )
                        # Add reference sequence
                        self._seq_append(
                            final_seq,
                            seqs[(i - 1) // 2][0][last_index : last_index + fill],
                            "R",
                            tuple(),
                            seqs[(i - 1) // 2][1][0] + last_index
                                + (fill - 1 if self.rev_strand else 0),
                            merge=False,
                        )
                        # Add edits
                        for edit in new_edits[pos_to_add]:
                            if edit[1] == "V" or edit[1] == "E":
                                # Edit is a substitution - set snv w/ variant info and insertion w/ placeholder
                                snv = (edit[0], edit[2], edit[3], edit[3][1])
                                insertion = (
                                    "",
                                    "R",
                                    tuple(),
                                    seqs[(i - 1) // 2][1][0] + fill,
                                )
                            else:
                                # Edit is an insertion
                                assert edit[1] == "I"
                                ins_index = bisect.bisect_left(
                                    [x[0] for x in intervals], edit[3][1] - 1
                                )
                                if (ins_index % 2) or (
                                    intervals[ins_index][0] == (edit[3][1] - 1)
                                ):
                                    # More reference sequence is needed to fill in - use snv for this
                                    snv = (
                                        self.bowtie_reference_index.get_stretch(
                                            self.chrom, edit[3][1] - 1, 1
                                        ),
                                        "R",
                                        tuple(),
                                        edit[3][1],
                                    )
                                else:
                                    # No reference sequence needed to fill in - use placeholder snv
                                    snv = ("", "R", tuple(), edit[3][1])
                                insertion = (edit[0], edit[2], edit[3], edit[3][1])
                        self._seq_append(final_seq, *snv, merge=False)
                        self._seq_append(final_seq, *insertion, merge=False)
                        last_index += fill + 1
                        last_pos += fill + 1
                    if intervals[i - 1][1] != "R":
                        # Exonic section bounded by deletion on left - add to seq
                        if isinstance(intervals[i - 1][2], list):
                            genomic_position = min([x[1] for x in intervals[i - 1][2]])
                        else:
                            genomic_position = intervals[i - 1][2][1]
                        if (
                            "",
                            intervals[i - 1][1],
                            [intervals[i - 1][2]],
                            genomic_position,
                        ) not in final_seq:
                            self._seq_append(
                                final_seq,
                                "",
                                intervals[i - 1][1],
                                intervals[i - 1][2],
                                genomic_position,
                                merge=False,
                            )
                    # Add reference sequence if necessary
                    ref_to_add = seqs[(i - 1) // 2][0][last_index:]
                    if ref_to_add:
                        if self.rev_strand:
                            self._seq_append(
                                final_seq,
                                ref_to_add,
                                "R",
                                tuple(),
                                seqs[(i - 1) // 2][1][1],
                                merge=False,
                            )
                        else:
                            self._seq_append(
                                final_seq,
                                ref_to_add,
                                "R",
                                tuple(),
                                seqs[(i - 1) // 2][1][0] + last_index,
                                merge=False,
                            )
                    if intervals[i][1] != "R":
                        # Exonic section bounded by deletion on right - add to seq
                        if isinstance(intervals[i][2], list):
                            genomic_position = min([x[1] for x in intervals[i][2]])
                        else:
                            genomic_position = intervals[i][2][1]
                        self._seq_append(
                            final_seq,
                            "",
                            intervals[i][1],
                            intervals[i][2],
                            genomic_position,
                            merge=False,
                        )
                    # Move to next exonic section
                    i += 2
                    try:
                        while pos > intervals[i][0]:
                            # Position is still outside of this exon
                            if intervals[i - 1][1] != "R":
                                # Exonic section bounded by deletion on left - add to seq
                                if isinstance(intervals[i - 1][2], list):
                                    genomic_position = min(
                                        [x[1] for x in intervals[i - 1][2]]
                                    )
                                else:
                                    genomic_position = intervals[i - 1][2][1]
                                self._seq_append(
                                    final_seq,
                                    "",
                                    intervals[i - 1][1],
                                    intervals[i - 1][2],
                                    genomic_position,
                                    merge=False,
                                )
                            # Add reference sequence
                            if self.rev_strand:
                                self._seq_append(
                                    final_seq,
                                    seqs[(i - 1) // 2][0],
                                    "R",
                                    tuple(),
                                    seqs[(i - 1) // 2][1][1],
                                    merge=False,
                                )
                            else:
                                self._seq_append(
                                    final_seq,
                                    seqs[(i - 1) // 2][0],
                                    "R",
                                    tuple(),
                                    seqs[(i - 1) // 2][1][0],
                                    merge=False,
                                )
                            if intervals[i][1] != "R":
                                # Exonic section bounded by deletion on right - add to seq
                                if isinstance(intervals[i][2], list):
                                    genomic_position = min(
                                        [x[1] for x in intervals[i][2]]
                                    )
                                else:
                                    genomic_position = intervals[i][2][1]
                                self._seq_append(
                                    final_seq,
                                    "",
                                    intervals[i][1],
                                    intervals[i][2],
                                    genomic_position,
                                    merge=False,
                                )
                            # Move to next exonic section
                            i += 2
                    except IndexError:
                        if i > len(intervals) - 1:
                            # Done enumerating sequence
                            break
                    pos_group = [pos]
                else:
                    pos_group.append(pos)
            if self.rev_strand:
                # Reverse complement sequence
                final_seq = [
                    (
                        seq[::-1].translate(revcomp_translation_table),
                        mutation_class,
                        orig_seq,
                        position,
                    )
                    for seq, mutation_class, orig_seq, position in final_seq
                ][::-1]
            adj_seq = [final_seq[0]]
            strand = 1 if not self.rev_strand else 0
            for i in range(1, len(final_seq)):
                if adj_seq[-1][0] == "" and final_seq[i][0] == "":
                    if len(adj_seq[-1][2]) == 2 and type(adj_seq[-1][2][0]) == list:
                        prev_pos = adj_seq[-1][3] + len(adj_seq[-1][2][0][2]) - 1
                    else:
                        prev_pos = (
                            adj_seq[-1][3]
                            + strand * sum([len(x[2]) for x in adj_seq[-1][2]])
                            - 1 * strand
                        )
                    if (
                        self.rev_strand
                        and (
                            final_seq[i][3]
                            + sum([len(x[2]) for x in final_seq[i][2]])
                            - 1
                        )
                        >= prev_pos
                    ) or (not self.rev_strand and final_seq[i][3] <= prev_pos):
                        adj_seq[-1] = self.hybridize_seq(
                            adj_seq[-1],
                            final_seq[i],
                            include_somatic=include_somatic,
                            include_germline=include_germline,
                        )
                    else:
                        adj_seq.append(final_seq[i])
                else:
                    adj_seq.append(final_seq[i])
            if not self.boundary_spanning_deletion:
                # Don't need to truncate the list of sequences
                for i in range(0, len(adj_seq)):
                    if adj_seq[i][0] == "":
                        if adj_seq[i][1] in ["G", "S"]:
                            adj_seq[i] = (
                                adj_seq[i][0],
                                adj_seq[i][1],
                                [adj_seq[i][2][0][:-1]],
                                adj_seq[i][3],
                            )
                        elif adj_seq[i][1] == "H":
                            adj_seq[i] = (
                                adj_seq[i][0],
                                adj_seq[i][1],
                                [
                                    adj_seq[i][2][0][0:4]
                                    + [[x[:-1] for x in adj_seq[i][2][0][4]]],
                                    adj_seq[i][2][1][0:4]
                                    + [[x[:-1] for x in adj_seq[i][2][1][4]]],
                                ],
                                adj_seq[i][3],
                            )
                return adj_seq
            else:
                # Need to determine where boundary spanning deletions start
                dels = [x for x in adj_seq if x[0] == ""]
                triggering_mutations = []
                for i in range(0, len(dels)):
                    if dels[i][1] in ["G", "S"]:
                        if dels[i][2][0][6]:
                            triggering_mutations.append(dels[i][2][0])
                            deletion_start = dels[i][3] - 1
                            start_index = bisect.bisect_left(
                                self.intervals, deletion_start
                            )
                            deletion_end = (dels[i][3] - 1) + len(dels[i][2][0][2]) - 1
                            end_index = bisect.bisect_left(self.intervals, deletion_end)
                            break
                    elif dels[i][1] == "H":
                        contributing_deletions = dels[i][2][0][4]
                        if len([x for x in contributing_deletions if x[6]]) > 0:
                            for x in contributing_deletions:
                                triggering_mutations.append(x)
                            deletion_start = dels[i][3] - 1
                            start_index = bisect.bisect_left(
                                self.intervals, deletion_start
                            )
                            deletion_end = (dels[i][3] - 1) + len(dels[i][2][0][2]) - 1
                            end_index = bisect.bisect_left(self.intervals, deletion_end)
                            break
                try:
                    if start_index % 2:
                        if self.rev_strand:
                            try:
                                acceptable = self.intervals[start_index + 1] + 2
                            except IndexError:
                                truncated_seq = []
                            else:
                                truncated_seq = [
                                    x
                                    for x in adj_seq
                                    if x[3] >= acceptable
                                    and x[3] <= self.intervals[-1] + 1
                                ]
                        else:
                            try:
                                acceptable = self.intervals[max(start_index - 2, 0)] + 1
                            except IndexError:
                                truncated_seq = []
                            else:
                                truncated_seq = [
                                    x
                                    for x in adj_seq
                                    if x[3] <= acceptable
                                    and x[3] >= self.intervals[0] + 2
                                ]
                    elif end_index % 2:
                        if self.rev_strand:
                            try:
                                acceptable = self.intervals[end_index + 1] + 2
                            except IndexError:
                                truncated_seq = []
                            else:
                                truncated_seq = [
                                    x
                                    for x in adj_seq
                                    if x[3] >= acceptable
                                    and x[3] <= self.intervals[-1] + 1
                                ]
                        else:
                            try:
                                acceptable = self.intervals[max(end_index - 2, 0)] + 1
                            except IndexError:
                                truncated_seq = []
                            else:
                                truncated_seq = [
                                    x
                                    for x in adj_seq
                                    if x[3] <= acceptable
                                    and x[3] >= self.intervals[0] + 2
                                ]
                    skipped_mutations = [
                        x for x in adj_seq if x[1] != "R" and x not in truncated_seq
                    ]
                    warnings.warn(
                        "".join(
                            [
                                "Deletion spanning intron-exon boundary",
                                " leads to truncated version of transcript ",
                                self.transcript_id,
                                ";\nContributing deletions: ",
                                str(triggering_mutations),
                                ";\nMutations excluded from trancript: ",
                                str(skipped_mutations),
                            ]
                        )
                    )
                except UnboundLocalError:
                    truncated_seq = [x for x in adj_seq]
                for i in range(0, len(truncated_seq)):
                    if truncated_seq[i][0] == "":
                        if truncated_seq[i][1] in ["G", "S"]:
                            truncated_seq[i] = (
                                truncated_seq[i][0],
                                truncated_seq[i][1],
                                [truncated_seq[i][2][0][:-1]],
                                truncated_seq[i][3],
                            )
                        elif truncated_seq[i][1] == "H":
                            truncated_seq[i] = (
                                truncated_seq[i][0],
                                truncated_seq[i][1],
                                [
                                    truncated_seq[i][2][0][0:4]
                                    + [[x[:-1] for x in truncated_seq[i][2][0][4]]],
                                    truncated_seq[i][2][1][0:4]
                                    + [[x[:-1] for x in truncated_seq[i][2][1][4]]],
                                ],
                                truncated_seq[i][3],
                            )
                return truncated_seq
        raise NotImplementedError(
            "Retrieving sequence with transcript coordinates not "
            "yet fully supported."
        )

    def _build_sequences(
        self, annotated_sequence, strand, include_somatic, include_germline, include_rna_edits
    ):
        """Builds alternative and reference sequences and dictionaries linking their
        transcript-level coordinates to genomic coordinates

        annotated_sequence: output of the annotated_seq() method, above
        strand: -1 for reverse strand transcript, 1 for forward strand

        Return value: alternative transcript sequence, reference transcript sequence,
            dictionary linking genomic coordinates (1-based) to reference
            transcript coordinates (0-based), dictionary linking genomic coordinates
            (1-based) to alternative transcript coordinates (0-based),
            dictionary linking reference transcript coordinates (0-based), to
            genomic coordinates (1-based), dictionary linking alternative transcript
            coordinates (0-based), to genomic coordinates (1-based)
        """
        counter, ref_counter = 0, 0  # hold edited transcript level coordinates
        genome_to_alt, genome_to_ref, alt_to_genome, ref_to_genome = (
            {},
            {},
            {},
            {},
        )  # hold coordinate linkers
        mut_to_ref_counter, mut_to_alt_counter = (
            {},
            {},
        )  # hold mut to transcript position info
        sequence, ref_sequence = "", ""  # hold flattened nucleotide sequence
        ref_tree, alt_tree = IntervalTree(), IntervalTree()
        # Process through sequence chunks to build reference/edited sequences
        # Link transcriptomic and genomic coordinates
        for seq in annotated_sequence:
            if seq[1] == "R":
                # Add sequence to both transcript versions
                sequence += seq[0]
                ref_sequence += seq[0]
                # Link coordinates
                genomic_positions = [seq[3] + (i * strand) for i in range(len(seq[0]))]
                transcriptomic_positions = [counter + i for i in range(len(seq[0]))]
                ref_transcriptomic_positions = [
                    ref_counter + i for i in range(len(seq[0]))
                ]
                multiassign(
                    genome_to_ref, genomic_positions, ref_transcriptomic_positions
                )
                multiassign(genome_to_alt, genomic_positions, transcriptomic_positions)
                multiassign(
                    ref_to_genome, ref_transcriptomic_positions, genomic_positions
                )
                multiassign(alt_to_genome, transcriptomic_positions, genomic_positions)
                # Update counters
                counter += len(seq[0])
                ref_counter += len(seq[0])
            elif seq[1] == "H":
                # Add relevant sequence to both transcripts
                sequence += seq[2][0][3]
                if self.rev_strand:
                    ref_sequence += seq[2][1][3][::-1].translate(
                        revcomp_translation_table
                    )
                else:
                    ref_sequence += seq[2][1][3]
                # Link ref and genomic coordinates
                deleted_length = len("".join([x[2] for x in seq[2][1][4]]))
                genomic_positions = [
                    (seq[3] + ((deleted_length - 1) * self.rev_strand)) + (i * strand)
                    for i in range(deleted_length)
                ]
                multiassign(
                    genome_to_ref,
                    genomic_positions,
                    [ref_counter for i in range(deleted_length)],
                )
                # Add variants to ref tree
                ref_tree[seq[3] : seq[3] + deleted_length] = seq[2][1][4]
                # Link alt and genomic coordinates
                deleted_length = len(seq[2][0][2])
                genomic_positions = [
                    (seq[3] + ((deleted_length - 1) * self.rev_strand)) + (i * strand)
                    for i in range(deleted_length)
                ]
                multiassign(
                    genome_to_alt,
                    genomic_positions,
                    [counter for i in range(deleted_length)],
                )
                # Add variants to alt tree
                alt_tree[seq[3] : seq[3] + deleted_length] = seq[2][0][4]
                # Add info to alt counter dict
                multiassign(
                    mut_to_alt_counter,
                    seq[2][0][4],
                    [(counter, seq) for i in range(len(seq[2][0][4]))],
                )
                # Update counters
                counter += len(seq[2][0][3])
                ref_counter += len(seq[2][1][3])
            elif seq[2][0][4] == "D":
                # Get length of deleted sequence
                deleted_length = len("".join([x[2] for x in seq[2]]))
                if not (
                    (seq[1] == "G" and include_germline == 2)
                    or (seq[1] == "S" and include_somatic == 2)
                ):
                    # Update reference transcript sequence + coordinates
                    variants = seq[2]
                    variants.sort(key=itemgetter(1), reverse=self.rev_strand)
                    for var in variants:
                        # Add sequence to reference transcript
                        if self.rev_strand:
                            ref_sequence += var[2][::-1].translate(
                                revcomp_translation_table
                            )
                        else:
                            ref_sequence += var[2]
                        # Link coordinates
                        # genomic_positions = [var[1] + (i*strand) for i in range(len(var[2]))]
                        genomic_positions = [
                            (var[1] + ((len(var[2]) - 1) * self.rev_strand))
                            + (i * strand)
                            for i in range(len(var[2]))
                        ]
                        transcriptomic_positions = [
                            ref_counter + i for i in range(len(var[2]))
                        ]
                        multiassign(
                            genome_to_ref, genomic_positions, transcriptomic_positions
                        )
                        multiassign(
                            ref_to_genome, transcriptomic_positions, genomic_positions
                        )
                        # Update counter
                        ref_counter += len(var[2])
                    # Update alternative coordinates
                    genomic_positions = [
                        (seq[3] + ((deleted_length - 1) * self.rev_strand))
                        + (i * strand)
                        for i in range(deleted_length)
                    ]
                    multiassign(
                        genome_to_alt,
                        genomic_positions,
                        [counter for i in range(deleted_length)],
                    )
                    # Update alt tree
                    alt_tree[seq[3] : seq[3] + deleted_length] = seq[2]
                    # Update alt counter info
                    multiassign(
                        mut_to_alt_counter,
                        seq[2],
                        [(counter, seq) for i in range(len(seq[2]))],
                    )
                else:
                    # Update ref coordinates only
                    genomic_positions = [
                        seq[3] + ((deleted_length - 1) * self.rev_strand) + (i * strand)
                        for i in range(deleted_length)
                    ]
                    multiassign(
                        genome_to_ref,
                        genomic_positions,
                        [ref_counter for i in range(deleted_length)],
                    )
                    multiassign(
                        genome_to_alt,
                        genomic_positions,
                        [counter for i in range(deleted_length)],
                    )
                    # Update ref tree
                    ref_tree[seq[3] : seq[3] + deleted_length] = seq[2]
            elif seq[2][0][4] == "I":
                # Add sequence to new transcript
                sequence += seq[0]
                # Link coordinates
                transcriptomic_positions = [counter + i for i in range(len(seq[0]))]
                multiassign(
                    alt_to_genome,
                    transcriptomic_positions,
                    [seq[3] for i in range(len(seq[0]))],
                )
                # Update alt tree
                alt_tree[seq[3] : seq[3] + 1] = seq[2]
                # Update alt counter info
                multiassign(
                    mut_to_alt_counter,
                    seq[2],
                    [(counter, seq) for i in range(len(seq[2]))],
                )
                # Update counter
                counter += len(seq[0])
                if (seq[1] == "G" and include_germline == 2) or (
                    seq[1] == "S" and include_somatic == 2
                ):
                    # Add sequence to reference transcript
                    ref_sequence += seq[0]
                    # Link coordinates
                    transcriptomic_positions = [
                        ref_counter + i for i in range(len(seq[0]))
                    ]
                    multiassign(
                        ref_to_genome,
                        transcriptomic_positions,
                        [seq[3] for i in range(len(seq[0]))],
                    )
                    # Update ref tree
                    ref_tree[seq[3] : seq[3] + 1] = seq[2]
                    # Update counter
                    ref_counter += len(seq[0])
            elif seq[2][0][4] == "V":
                # Add sequence to transcripts
                sequence += seq[0]
                if (seq[1] == "G" and include_germline == 2) or (
                    seq[1] == "S" and include_somatic == 2
                ):
                    ref_sequence += seq[0]
                    if self.rev_strand:
                        ref_tree[seq[3] - len(seq[0]) : seq[3]] = seq[2]
                    else:
                        ref_tree[seq[3] : seq[3] + len(seq[0])] = seq[2]
                else:
                    if self.rev_strand:
                        for i in seq[2]:
                            ref_sequence += i[2][::-1].translate(
                                revcomp_translation_table
                            )
                    else:
                        for i in seq[2]:
                            ref_sequence += i[2]
                # Link coordinates
                genomic_positions = [seq[3] + (i * strand) for i in range(len(seq[0]))]
                transcriptomic_positions = [counter + i for i in range(len(seq[0]))]
                ref_transcriptomic_positions = [
                    ref_counter + i for i in range(len(seq[0]))
                ]
                multiassign(
                    genome_to_ref, genomic_positions, ref_transcriptomic_positions
                )
                multiassign(genome_to_alt, genomic_positions, transcriptomic_positions)
                multiassign(
                    ref_to_genome, ref_transcriptomic_positions, genomic_positions
                )
                multiassign(alt_to_genome, transcriptomic_positions, genomic_positions)
                # Update alt trees
                if self.rev_strand:
                    alt_tree[seq[3] - len(seq[0]) : seq[3]] = seq[2]
                else:
                    alt_tree[seq[3] : seq[3] + len(seq[0])] = seq[2]
                # Update alt counter info
                multiassign(
                    mut_to_alt_counter,
                    seq[2],
                    [(counter, seq) for i in range(len(seq[2]))],
                )
                multiassign(
                    mut_to_ref_counter,
                    seq[2],
                    [ref_counter for i in range(len(seq[2]))],
                )
                # Update counters
                counter += len(seq[0])
                ref_counter += len(seq[0])
            elif seq[2][0][4] == "E":
                # Add RNA edits to transcripts
                sequence += seq[0]
                if include_rna_edits == 2:
                    ref_sequence += seq[0]
                    if self.rev_strand:
                        ref_tree[seq[3] - len(seq[0]): seq[3]] = seq[2]
                    else:
                        ref_tree[seq[3]: seq[3] + len(seq[0])] = seq[2]
                else:
                    if self.rev_strand:
                        for i in seq[2]:
                            ref_sequence += i[2][::-1].translate(
                                revcomp_translation_table
                            )
                    else:
                        for i in seq[2]:
                            ref_sequence += i[2]
                # Link coordinates
                genomic_positions = [seq[3] + (i * strand) for i in range(len(seq[0]))]
                transcriptomic_positions = [counter + i for i in range(len(seq[0]))]
                ref_transcriptomic_positions = [
                    ref_counter + i for i in range(len(seq[0]))
                ]
                multiassign(
                    genome_to_ref, genomic_positions, ref_transcriptomic_positions
                )
                multiassign(genome_to_alt, genomic_positions, transcriptomic_positions)
                multiassign(
                    ref_to_genome, ref_transcriptomic_positions, genomic_positions
                )
                multiassign(alt_to_genome, transcriptomic_positions, genomic_positions)
                # Update alt trees
                if self.rev_strand:
                    alt_tree[seq[3] - len(seq[0]) : seq[3]] = seq[2]
                else:
                    alt_tree[seq[3] : seq[3] + len(seq[0])] = seq[2]
                # Update alt counter info
                multiassign(
                    mut_to_alt_counter,
                    seq[2],
                    [(counter, seq) for i in range(len(seq[2]))],
                )
                multiassign(
                    mut_to_ref_counter,
                    seq[2],
                    [ref_counter for i in range(len(seq[2]))],
                )
                # Update counters
                counter += len(seq[0])
                ref_counter += len(seq[0])
        return (
            sequence,
            ref_sequence,
            genome_to_ref,
            genome_to_alt,
            ref_to_genome,
            alt_to_genome,
            ref_tree,
            alt_tree,
            mut_to_ref_counter,
            mut_to_alt_counter,
        )

    def neopeptides(
        self,
        min_size=8,
        max_size=11,
        include_somatic=1,
        include_germline=2,
        include_rna_edits=0,
        only_novel_upstream=False,
        only_downstream=True,
        only_reference=False,
        return_protein=False,
        allow_no_edits=False,
        allow_partial_codons=False
    ):
        """Retrieves dict of predicted peptide fragments from transcript that
        arise from one or more variants.
        min_size: minimum subpeptide length (specified as # of amino acids)
        max_size: maximum subpeptide length (specified as # of amino acids)
        include_somatic: 0 = do not include somatic mutations,
            1 = exclude somatic mutations from reference comparison,
            2 = include somatic mutations in both annotated sequence and
            reference comparison
        include_germline: 0 = do not include germline mutations,
            1 = exclude germline mutations from reference comparison,
            2 = include germline mutations in both annotated sequence and
            reference comparison
        include_rna_edits: 0 = do not include A to I RNA editing,
            1 = exclude RNA edits from reference comparison,
            2 = include RNA edits in both annotated sequence and reference comparison
        Return value: dict of peptides of desired length(s) [KEYS] with
            values equivalent to a list of causal variants [VALUES].
        """
        # if no edits to process, then skip all next steps and return {}
        if include_somatic == include_germline and include_somatic != 1:
            if not return_protein:
                return {}
            else:
                return {}, ""
        # Check whether to return due to lack of edits
        if not self.edits and not self.deletion_intervals and not allow_no_edits:
            if not return_protein:
                return {}
            else:
                return {}, ""
        # min size to process is 2 amino acids, otherwise skip and return {}
        if min_size < 2:
            if not return_protein:
                return {}
            else:
                return {}, ""
        # ensure max_size is not smaller than min_size
        if max_size < min_size:
            max_size = min_size
        # Pull nucleotide sequence
        annotated_seq = self.annotated_seq(
            include_somatic=include_somatic, include_germline=include_germline,
            include_rna_edits=include_rna_edits
        )
        # Check whether transcript is missing start codon
        if self.start_codon is None or "".join([seq[0] for seq in annotated_seq]) == "":
            if not return_protein:
                return {}
            else:
                return {}, ""
        # Set transcript variables
        transcript_warnings = []
        unknown_aa = False
        strand = 1 - self.rev_strand * 2  # +1 is + strand, -1 is - strand
        # Build alt/ref transcript sequences and generate position dictionaries/trees
        (
            sequence,
            ref_sequence,
            genome_to_ref,
            genome_to_alt,
            ref_to_genome,
            alt_to_genome,
            ref_tree,
            alt_tree,
            mut_to_ref_counter,
            mut_to_alt_counter,
        ) = self._build_sequences(
            annotated_seq,
            strand=strand,
            include_somatic=include_somatic,
            include_germline=include_germline,
            include_rna_edits=include_rna_edits
        )
        # Assign possible start codon sequences
        # May want to update to include 'ATA', 'ATT' for mitochondrial transcripts
        start_seqs = [self.start_codon_seq, "ATG"]
        # Check ref transcript start codon sequence
        ref_start_disrupted = False
        try:
            ref_tx_start = genome_to_ref[self.start_codon] - 2 * self.rev_strand
        except KeyError:
            if not only_reference:
                ref_start_disrupted = True
                warnings.warn(
                    "".join(
                        [
                            "Annotated start codon disrupted for transcript ",
                            self.transcript_id,
                            " by background mutation",
                        ]
                    ),
                    Warning,
                )
                transcript_warnings.append("annotated_start_codon_disrupted_background")
                coords = list(genome_to_ref.keys())
                coords.sort(reverse=self.rev_strand)
                start_index = min(
                    bisect.bisect_left(coords, self.start_codon), len(coords) - 1
                )
                ref_tx_start = coords[start_index]
            else:
                if not return_protein:
                    return {}
                else:
                    return {}, ""
        else:
            if ref_sequence[ref_tx_start : ref_tx_start + 3] not in start_seqs:
                warnings.warn(
                    "".join(
                        [
                            "Annotated start codon disrupted for transcript ",
                            self.transcript_id,
                            " by background mutation",
                        ]
                    ),
                    Warning,
                )
                # A mutation disrupted the annotated start codon
                if only_reference:
                    # Not compatible with user settings
                    if not return_protein:
                        return {}
                    else:
                        return {}, ""
                ref_start_disrupted = True
                transcript_warnings.append("annotated_start_codon_disrupted_background")
        # Check alt transcript start codon sequence
        alt_start_disrupted = False
        try:
            alt_tx_start = genome_to_alt[self.start_codon] - 2 * self.rev_strand
        except KeyError:
            if not only_reference:
                alt_start_disrupted = True
                warnings.warn(
                    "".join(
                        [
                            "Annotated start codon disrupted for transcript ",
                            self.transcript_id,
                        ]
                    ),
                    Warning,
                )
                transcript_warnings.append("annotated_start_codon_disrupted")
                coords = list(genome_to_alt.keys())
                coords.sort(reverse=self.rev_strand)
                start_index = min(
                    bisect.bisect_left(coords, self.start_codon), len(coords) - 1
                )
                alt_tx_start = coords[start_index]
            else:
                if not return_protein:
                    return {}
                else:
                    return {}, ""
        else:
            if sequence[alt_tx_start : alt_tx_start + 3] not in start_seqs:
                # A mutation disrupted the annotated start codon
                if not ref_start_disrupted:
                    warnings.warn(
                        "".join(
                            [
                                "Annotated start codon disrupted for transcript ",
                                self.transcript_id,
                            ]
                        ),
                        Warning,
                    )
                    transcript_warnings.append("annotated_start_codon_disrupted")
                if only_reference:
                    # Not compatible with user settings
                    if not return_protein:
                        return {}
                    else:
                        return {}, ""
                alt_start_disrupted = True
        # Set reference for atg sequences
        ref_atg = None
        atg_positions = [m.start() for m in re.finditer("ATG", ref_sequence)]
        # First check for upstream ATGs if relevant
        if not only_reference and not only_downstream:
            # find start index
            index = bisect.bisect_left(atg_positions, ref_tx_start)
            for atg in atg_positions[0:index]:
                # Start codon is before annotated start - get reference genome seq
                genome_pos = [
                    ref_to_genome[atg],
                    ref_to_genome[atg + 1],
                    ref_to_genome[atg + 2],
                ]
                ref_genome_seq = "".join(
                    [
                        self.bowtie_reference_index.get_stretch(self.chrom, x - 1, 1)
                        for x in genome_pos
                    ]
                )
                if self.rev_strand:
                    ref_genome_seq = ref_genome_seq[::-1].translate(
                        revcomp_translation_table
                    )
                if not only_novel_upstream or (
                    only_novel_upstream
                    and ref_sequence[atg : atg + 3] != ref_genome_seq
                ):
                    # Valid start codon for user settings
                    alt_equiv = genome_to_alt[genome_pos[0]]
                    muts = set()
                    for pos in genome_pos:
                        muts.update([x.data[0] for x in ref_tree[pos : pos + 1]])
                    ref_atg = [alt_equiv, atg, list(muts)]
                    break
        # Check for annotated start if not using upstream start
        if not ref_start_disrupted and ref_atg is None:
            # Annotated start undisrupted
            alt_equiv = genome_to_alt[ref_to_genome[ref_tx_start]]
            ref_atg = [alt_equiv, ref_tx_start, []]
        elif ref_atg is None:
            # Get mutations overlapping start position
            muts = set()
            for pos in self.start_coordinates:
                muts.update([x.data[0] for x in ref_tree[pos : pos + 1]])
            # Find new start if available
            index = bisect.bisect_left(atg_positions, ref_tx_start)  # find start index
            for atg in atg_positions[index:]:
                genome_pos = [
                    ref_to_genome[atg],
                    ref_to_genome[atg + 1],
                    ref_to_genome[atg + 2],
                ]
                # Downstream start
                alt_equiv = genome_to_alt[genome_pos[0]]
                for pos in genome_pos:
                    muts.update([x.data[0] for x in ref_tree[pos : pos + 1]])
                ref_atg = [alt_equiv, atg, list(muts)]
                break
        # Set alternative atg
        start_codon = None
        atg_positions = [m.start() for m in re.finditer("ATG", sequence)]
        # First check for upstream ATGs if relevant
        if not only_downstream and not only_reference:
            # find start index
            index = bisect.bisect_left(atg_positions, alt_tx_start)
            for atg in atg_positions[0:index]:
                # Start codon is before annotated start - get reference tx seq
                genome_pos = [
                    alt_to_genome[atg],
                    alt_to_genome[atg + 1],
                    alt_to_genome[atg + 2],
                ]
                ref_tx_pos = genome_to_ref[genome_pos[0]]
                ref_tx_seq = ref_sequence[ref_tx_pos : ref_tx_pos + 3]
                if not only_novel_upstream or (
                    only_novel_upstream and sequence[atg : atg + 3] != ref_tx_seq
                ):
                    # Valid start codon for user settings
                    muts = set()
                    for pos in genome_pos:
                        muts.update([x.data[0] for x in alt_tree[pos : pos + 1]])
                    start_codon = [atg, ref_tx_pos, list(muts)]
                    # Assess reading frame offset to reference tx
                    try:
                        reading_frame = (start_codon[0] - ref_atg[0]) % 3
                    except TypeError:
                        reading_frame = None
                    # Give warning
                    if sequence[atg : atg + 3] != ref_tx_seq:
                        novelty = "novel"
                    else:
                        novelty = "preexisting"
                    warnings.warn(
                        " ".join(
                            [
                                "Using",
                                novelty,
                                "start codon",
                                "upstream of",
                                "reference start codon for",
                                self.transcript_id,
                            ]
                        )
                    )
                    transcript_warnings.append(
                        "_".join(
                            [
                                novelty,
                                "start",
                                "codon",
                                "upstream",
                                "of",
                                "reference",
                                "start",
                                "codon",
                            ]
                        )
                    )
                    break
        # Check for annotated start if not using upstream start
        if not alt_start_disrupted and start_codon is None:
            ref_equiv = genome_to_ref[alt_to_genome[alt_tx_start]]
            start_codon = [alt_tx_start, ref_equiv, []]
            try:
                reading_frame = (start_codon[0] - ref_atg[0]) % 3
            except TypeError:
                reading_frame = None
        elif start_codon is None:
            # Get mutations overlapping start position
            muts = set()
            for pos in self.start_coordinates:
                muts.update([x.data[0] for x in alt_tree[pos : pos + 1]])
            # Find new start if available
            index = bisect.bisect_left(atg_positions, alt_tx_start)  # find start index
            for atg in atg_positions[index:]:
                genome_pos = [
                    alt_to_genome[atg],
                    alt_to_genome[atg + 1],
                    alt_to_genome[atg + 2],
                ]
                if atg > alt_tx_start:
                    # Downstream start
                    ref_equiv = genome_to_ref[genome_pos[0]]
                    for pos in genome_pos:
                        muts.update([x.data[0] for x in alt_tree[pos : pos + 1]])
                    start_codon = [atg, ref_equiv, list(muts)]
                    # Assess reading frame offset to reference tx
                    try:
                        reading_frame = (start_codon[0] - ref_atg[0]) % 3
                    except TypeError:
                        reading_frame = None
                    # Give warning
                    if (
                        sequence[atg : atg + 3]
                        != ref_sequence[ref_equiv : ref_equiv + 3]
                    ):
                        novelty = "novel"
                    else:
                        novelty = "preexisting"
                    warnings.warn(
                        " ".join(
                            [
                                "Using",
                                novelty,
                                "start codon",
                                "downstream of",
                                "reference start codon for",
                                self.transcript_id,
                            ]
                        )
                    )
                    transcript_warnings.append(
                        "_".join(
                            [
                                novelty,
                                "start",
                                "codon",
                                "downstream",
                                "of",
                                "reference",
                                "start",
                                "codon",
                            ]
                        )
                    )
                    break
        # Check whether we have valid start codon for alt sequence - if not, return nothing
        if start_codon is None:
            warnings.warn(
                "".join(
                    [
                        "Start codon disrupted for transcript ",
                        self.transcript_id,
                        "; no valid peptides",
                    ]
                ),
                Warning,
            )
            if not return_protein:
                return {}
            else:
                return {}, ""
        # Check whether using same start codon (or at least frame) for ref/alternative transcripts
        same_start_frame = True
        if ref_atg is not None:
            # Translation occurs on ref transcript
            if alt_to_genome[start_codon[0]] != ref_to_genome[ref_atg[1]]:
                # Different genomic positions for start codons
                if (start_codon[0] - ref_atg[0]) % 3:
                    # Different reading frames between the transcripts
                    same_start_frame = False
        else:
            # No translation from ref transcript - "new" reading frame
            same_start_frame = False
        # Assign valid stop codon sequences
        if not self.mitochondrial:
            stop_seqs = [self.stop_codon_seq, "TAA", "TGA", "TAG"]
        else:
            stop_seqs = [self.stop_codon_seq, "TAA", "TAG", "AGA", "AGG"]
        ref_stop_disrupted = False
        alt_stop_disrupted = False
        if self.stop_codon is not None:
            # Check ref transcript stop codon sequence if translation occurs
            if ref_atg is not None:
                try:
                    ref_tx_stop = genome_to_ref[self.stop_codon] - 2 * self.rev_strand
                except KeyError:
                    ref_stop_disrupted = True
                    warnings.warn(
                        "".join(
                            [
                                "Annotated stop codon disrupted for transcript ",
                                self.transcript_id,
                                " by background mutation",
                            ]
                        ),
                        Warning,
                    )
                    transcript_warnings.append(
                        "annotated_stop_codon_disrupted_background"
                    )
                else:
                    if ref_sequence[ref_tx_stop : ref_tx_stop + 3] not in stop_seqs:
                        # Mutation disrupted stop codon
                        ref_stop_disrupted = True
                    elif ref_atg[1] > ref_tx_stop:
                        # Translation starts downstream of annotated stop codon
                        ref_stop_disrupted = True
                    elif (ref_tx_stop - ref_atg[1]) % 3:
                        # Annotated stop is out of frame of reference tx start
                        ref_stop_disrupted = True
                    if ref_stop_disrupted and not (
                        not ref_start_disrupted and (ref_tx_stop - ref_atg[1]) % 3
                    ):
                        warnings.warn(
                            "".join(
                                [
                                    "Annotated stop codon disrupted for transcript ",
                                    self.transcript_id,
                                    " by background mutation",
                                ]
                            ),
                            Warning,
                        )
                        transcript_warnings.append(
                            "annotated_stop_codon_disrupted_background"
                        )
            # Check alt transcript start codon sequence
            try:
                alt_tx_stop = genome_to_alt[self.stop_codon] - 2 * self.rev_strand
            except KeyError:
                alt_stop_disrupted = True
                transcript_warnings.append("annotated_stop_codon_disrupted")
                warnings.warn(
                    "".join(
                        [
                            "Annotated stop codon disrupted for transcript ",
                            self.transcript_id,
                        ]
                    ),
                    Warning,
                )
            else:
                if sequence[alt_tx_stop : alt_tx_stop + 3] not in stop_seqs:
                    # Mutation disrupted stop codon
                    alt_stop_disrupted = True
                elif start_codon[0] > alt_tx_stop:
                    # Translation starts downstream of annotated stop codon
                    alt_stop_disrupted = True
                elif (alt_tx_stop - start_codon[0]) % 3:
                    # Annotated stop is out of frame of alt tx start
                    alt_stop_disrupted = True
                if (
                    alt_stop_disrupted
                    and not ref_stop_disrupted
                    and not (
                        not alt_start_disrupted and (alt_tx_stop - start_codon[0]) % 3
                    )
                ):
                    warnings.warn(
                        "".join(
                            [
                                "Annotated stop codon disrupted for transcript ",
                                self.transcript_id,
                            ]
                        ),
                        Warning,
                    )
                    transcript_warnings.append("annotated_stop_codon_disrupted")
        else:
            transcript_warnings.append("no_annotated_stop_codon")
        # Set up reference stop
        ref_stop = None
        if ref_atg is not None:
            muts = set()
            if self.stop_codon is not None:
                # Check for reference tx stop codon if translation occurs
                if not ref_stop_disrupted:
                    # Annotated stop codon could be valid - not disrupted and in-frame
                    alt_equiv = genome_to_alt[ref_to_genome[ref_tx_stop]]
                    ref_stop = [alt_equiv, ref_tx_stop, []]
                else:
                    # Find variants disrupting stop codon if applicable
                    for pos in self.stop_coordinates:
                        muts.update([x.data[0] for x in ref_tree[pos : pos + 1]])
            # Find new stop codon if relevant
            stop_positions = []
            for seq in stop_seqs[1:]:
                # Save all possible stop positions that are in frame
                stop_positions.extend(
                    [
                        m.start()
                        for m in re.finditer(seq, ref_sequence)
                        if m.start() > (ref_atg[1] + 2)
                        and not ((m.start() - ref_atg[1]) % 3)
                    ]
                )
            stop_positions.sort()
            annotated_start = self.start_coordinates[0 + 2 * self.rev_strand]
            indels = False
            if [
                x
                for x in annotated_seq
                if x[1] != "R" and (x[0] == "" or x[2][0][4] == "I")
            ]:
                indels = True
            for stop in stop_positions:
                genome_pos = [
                    ref_to_genome[stop],
                    ref_to_genome[stop + 1],
                    ref_to_genome[stop + 2],
                ]
                # Check novelty
                ref_genome_seq = "".join(
                    [
                        self.bowtie_reference_index.get_stretch(self.chrom, x - 1, 1)
                        for x in genome_pos
                    ]
                )
                if self.rev_strand:
                    ref_genome_seq = ref_genome_seq[::-1].translate(
                        revcomp_translation_table
                    )
                novel = False
                if ref_sequence[stop : stop + 3] != ref_genome_seq:
                    # Novel sequence via mutation
                    novel = True
                elif ref_start_disrupted or indels:
                    # Reading frame may be different - compare to annotated sequence
                    if bisect.bisect_left(
                        self.intervals, genome_pos[0]
                    ) == bisect.bisect_left(self.intervals, annotated_start):
                        # potential stop is in same exon as annotated start codon
                        if (genome_pos[0] - annotated_start) % 3:
                            novel = True
                    else:
                        # potential stop is in different exon
                        stop_index = bisect.bisect_left(
                            self.intervals, genome_pos[0] - 1
                        )
                        start_index = bisect.bisect_left(
                            self.intervals, annotated_start - 1
                        )
                        if stop_index > start_index:
                            seq_length = (
                                self.intervals[start_index] - annotated_start + 2
                            )
                            seq_length += (
                                genome_pos[0] - self.intervals[stop_index - 1] - 2
                            )
                            for i in range(start_index + 1, stop_index - 1, 2):
                                seq_length += self.intervals[i + 1] - self.intervals[i]
                        else:
                            seq_length = self.intervals[stop_index] - genome_pos[0] + 2
                            seq_length += (
                                annotated_start - self.intervals[start_index - 1] - 2
                            )
                            for i in range(stop_index + 1, start_index - 1, 2):
                                seq_length += self.intervals[i + 1] - self.intervals[i]
                        if seq_length % 3:
                            # Frame shift relative to annotated transcript
                            novel = True
                # Check whether to use stop
                use_stop = False
                if ref_stop is not None and stop < ref_stop[1] and novel:
                    # Stop codon is before intact reference stop and novel
                    use_stop = True
                elif ref_stop is None and self.stop_codon is not None:
                    # Reference stop is disrupted
                    if self.rev_strand and (
                        (genome_pos[0] < self.stop_codon + 2)
                        or (genome_pos[0] > self.stop_codon + 2 and novel)
                    ):
                        # Reverse strand - stop is after annotated stop or before and novel
                        use_stop = True
                    elif not self.rev_strand and (
                        (genome_pos[0] > self.stop_codon)
                        or (genome_pos[0] < self.stop_codon and novel)
                    ):
                        # Forward strand - stop is after annotated stop or before and novel
                        use_stop = True
                elif ref_stop is None and self.stop_codon is None and novel:
                    # No reference stop codon, this stop is novel
                    use_stop = True
                if use_stop:
                    # This is the stop codon to use
                    alt_equiv = genome_to_alt[genome_pos[0]]
                    for pos in genome_pos:
                        muts.update([x.data[0] for x in ref_tree[pos : pos + 1]])
                    ref_stop = [alt_equiv, stop, list(muts)]
                    break
        # Set up alternative stop
        stop_codon = None
        muts = set()
        if (
            ref_stop is not None
            and (not (ref_stop[0] - start_codon[0]) % 3)
            and sequence[ref_stop[0] : ref_stop[0] + 3] in stop_seqs
        ):
            stop_codon = [ref_stop[0], ref_stop[1], []]
        elif self.stop_codon is not None:
            if not alt_stop_disrupted:
                # Annotated stop codon could be valid - not disrupted and in-frame
                ref_equiv = genome_to_ref[alt_to_genome[alt_tx_stop]]
                stop_codon = [alt_tx_stop, ref_equiv, []]
            else:
                # Find variants disrupting stop codon if applicable
                for pos in self.stop_coordinates:
                    muts.update([x.data[0] for x in alt_tree[pos : pos + 1]])
        # Find new stop codon if relevant
        stop_positions = []
        for seq in stop_seqs[1:]:
            # save all possible stop positions that are in frame
            stop_positions.extend(
                [
                    m.start()
                    for m in re.finditer(seq, sequence)
                    if m.start() > (start_codon[0] + 2)
                    and not ((m.start() - start_codon[0]) % 3)
                ]
            )
        stop_positions.sort()
        for stop in stop_positions:
            genome_pos = [
                alt_to_genome[stop],
                alt_to_genome[stop + 1],
                alt_to_genome[stop + 2],
            ]
            # Check novelty
            ref_tx_seq = ref_sequence[
                genome_to_ref[genome_pos[0]] : genome_to_ref[genome_pos[0]] + 3
            ]
            novel = False
            if sequence[stop : stop + 3] != ref_tx_seq or (
                ref_atg is not None
                and ((genome_to_ref[genome_pos[0]] - ref_atg[1]) % 3)
            ):
                novel = True
            # Check whether to use stop
            use_stop = False
            if stop_codon is not None and stop < stop_codon[0] and novel:
                # Stop codon is before intact reference stop and novel
                use_stop = True
            elif stop_codon is None and self.stop_codon is not None:
                # Reference stop is disrupted
                if self.rev_strand and (
                    (genome_pos[0] < self.stop_codon + 2)
                    or (genome_pos[0] > self.stop_codon + 2 and novel)
                ):
                    # Reverse strand - stop is after annotated stop or before and novel
                    use_stop = True
                elif not self.rev_strand and (
                    (genome_pos[0] > self.stop_codon)
                    or (genome_pos[0] < self.stop_codon and novel)
                ):
                    # Forward strand - stop is after annotated stop or before and novel
                    use_stop = True
            elif stop_codon is None and self.stop_codon is None and novel:
                # No reference stop codon, this stop is novel
                use_stop = True
            if use_stop:
                # This is the stop codon to use
                ref_equiv = genome_to_ref[genome_pos[0]]
                for pos in genome_pos:
                    muts.update([x.data[0] for x in alt_tree[pos : pos + 1]])
                stop_codon = [stop, ref_equiv, list(muts)]
                break
        # Create shortcuts to translation start positions
        ref_start = start_codon[1]
        coding_start = start_codon[0]
        # Store coordinates of variants and frame shifts
        coordinates = []
        # Frame shifts: [genomic start coordinate, genomic end coordinate, CDS-
        # level start coordinate, CDS-level end coordinate, mutation info
        # associated with frame shift]
        frame_shifts = []
        if reading_frame:
            # New start codon is in different frame that ref
            frame_shifts.append(
                [alt_to_genome[coding_start], -1, 0, -1, start_codon[2]]
            )
        if alt_stop_disrupted and not reading_frame:
            # Alternative stop disrupted due to reasons besides new reading frame (variant)
            if ref_stop is not None and stop_codon is not None:
                # Both transcripts have valid stop codon
                if ref_stop[0] != stop_codon[0]:
                    # Occur at different genomic positions
                    frame_shifts.append(
                        [
                            alt_to_genome[ref_stop[0]],
                            alt_to_genome[stop_codon[0]],
                            ref_stop[0],
                            stop_codon[0],
                            stop_codon[2],
                        ]
                    )
            elif stop_codon is not None:
                # No reference stop but alt stop - adjust accordingly
                frame_shifts.append(
                    [
                        alt_to_genome[alt_tx_stop],
                        alt_to_genome[stop_codon[0]],
                        alt_tx_stop,
                        stop_codon[0],
                        stop_codon[2],
                    ]
                )
            else:
                # No valid stop in new transcripts - flag transcript with warning
                warnings.warn(
                    "".join(
                        [
                            "Stop codon not detected prior",
                            " to end of transcript ",
                            self.transcript_id,
                            "; this",
                            " transcript may undergo ",
                            "degradation",
                        ]
                    ),
                    Warning,
                )
                transcript_warnings.append("nonstop")
        elif stop_codon is not None and ref_stop is not None:
            # Both transcripts have stop and no disruptions
            if stop_codon[0] > ref_stop[0] and same_start_frame:
                # Stop codons in same frame, but new one occurs later than ref
                frame_shifts.append(
                    [
                        alt_to_genome[ref_stop[0]],
                        alt_to_genome[stop_codon[0]],
                        ref_stop[0],
                        stop_codon[0],
                        ref_stop[2],
                    ]
                )
        # Get coordinates for epitope enumeration
        if self.rev_strand:
            try:
                relevant_variants = [
                    x.data[0]
                    for x in alt_tree[
                        alt_to_genome[stop_codon[0]] : alt_to_genome[start_codon[0]] + 1
                    ]
                ]
            except TypeError:
                relevant_variants = [
                    x.data[0]
                    for x in alt_tree[
                        alt_to_genome[len(sequence) - 1] : alt_to_genome[start_codon[0]]
                        + 1
                    ]
                ]
        else:
            try:
                relevant_variants = [
                    x.data[0]
                    for x in alt_tree[
                        alt_to_genome[start_codon[0]] : alt_to_genome[stop_codon[0]]
                    ]
                ]
            except TypeError:
                relevant_variants = [
                    x.data[0]
                    for x in alt_tree[
                        alt_to_genome[start_codon[0]] : alt_to_genome[len(sequence) - 1]
                        + 1
                    ]
                ]
        for var in relevant_variants:
            # Get alt transcript counter position at variant and full annotated seq chunk
            counter, seq = mut_to_alt_counter[var]
            try:
                # Get ref counter at variant
                ref_counter = mut_to_ref_counter[var]
            except KeyError:
                # Not relevant for variant
                ref_counter = None
            if counter < coding_start + 2:
                if seq[1] == "H":
                    # Hybrid deletion overlapping start
                    if len(seq[2][0][3]) + counter >= coding_start:
                        coordinates.append(
                            [
                                alt_to_genome[coding_start],
                                seq[2][0][1] + len(seq[2][0][3]) * strand - 1,
                                counter,
                                counter + len(seq[2][0][3]) - 1,
                                "NA",
                                "NA",
                                seq[2][0][4],
                            ]
                        )
                    continue
                elif seq[2][0][4] == "D":
                    # Deletion overlapping start
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                alt_to_genome[coding_start],
                                seq[3] + len(seq[0]) * strand - 1,
                                counter,
                                counter + len(seq[0]) - 1,
                                "NA",
                                "NA",
                                seq[2],
                            ]
                        )
                    continue
                elif seq[2][0][4] == "I":
                    # Insertion overlapping start
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                alt_to_genome[coding_start],
                                seq[3] + len(seq[0]) * strand - 1,
                                0,
                                counter + len(seq[0]) - coding_start - 1,
                                "NA",
                                "NA",
                                seq[2],
                            ]
                        )
                    continue
                elif seq[2][0][4] == "V":
                    # SNV overlapping start
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                alt_to_genome[coding_start],
                                seq[3] + len(seq[0]) * strand - 1,
                                0,
                                counter + len(seq[0]) - coding_start - 1,
                                0,
                                ref_counter + len(seq[0]) - ref_start - 1,
                                seq[2],
                            ]
                        )
                    continue
                elif seq[2][0][4] == "E":
                    # RNA edit overlapping start
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                alt_to_genome[coding_start],
                                seq[3] + len(seq[0]) * strand - 1,
                                0,
                                counter + len(seq[0]) - coding_start - 1,
                                0,
                                ref_counter + len(seq[0]) - ref_start - 1,
                                seq[2],
                            ]
                        )
                    continue
            if seq[1] == "H":
                # Hybrid deletion downstream of start
                coordinates.append(
                    [
                        seq[3],
                        seq[3] + len(seq[0]) * strand - 1,
                        counter,
                        counter + len(seq[0]) - 1,
                        "NA",
                        "NA",
                        seq[2][0][4],
                    ]
                )
                read_frame2 = self.reading_frame(seq[2][0][1] + len(seq[2][0][3]))
                if read_frame2 is None:
                    # this case NOT addressed at present
                    # (e.g. deletion involves all or part of intron)
                    continue
                read_frame1 = (read_frame2 + len(seq[2][1][3])) % 3
                if read_frame1 != read_frame2:
                    # splicing variation (e.g. deletion of part of intron/exon)
                    if reading_frame == 0:
                        reading_frame = (read_frame1 - read_frame2) % 3
                        frame_shifts.append(
                            [seq[2][0][1], -1, counter, -1, seq[2][0][4]]
                        )
                    elif (reading_frame + read_frame1 - read_frame2) % 3 == 0:
                        # close out all frame_shifts ending in -1
                        for i in range(len(frame_shifts), 0, -1):
                            if frame_shifts[i - 1][1] < 0:
                                frame_shifts[i - 1][1] = seq[3] + len(seq[2][0][2])
                                frame_shifts[i - 1][3] = counter + len(seq[2][0][2])
                            else:
                                continue
                        reading_frame = 0
                    else:
                        frame_shifts.append(
                            [seq[2][0][1], -1, counter, -1, seq[2][0][4]]
                        )
                        reading_frame = (reading_frame + read_frame1 - read_frame2) % 3
                continue
            elif seq[2][0][4] == "D":
                # Deletion downstream of start
                coordinates.append(
                    [
                        seq[3],
                        seq[3] + len(seq[0]) * strand - 1,
                        counter,
                        counter + len(seq[0]) - 1,
                        "NA",
                        "NA",
                        seq[2],
                    ]
                )
                read_frame1 = self.reading_frame(seq[3])
                read_frame2 = self.reading_frame(seq[3] + len(seq[2][0][2]))
                if read_frame1 is None or read_frame2 is None:
                    # these cases NOT addressed at present
                    # (e.g. deletion involves all or part of intron)
                    continue
                if read_frame1 != read_frame2:
                    # splicing variation (e.g. deletion of part of intron/exon)
                    if reading_frame == 0:
                        reading_frame = (read_frame1 - read_frame2) % 3
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                    elif (reading_frame + read_frame1 - read_frame2) % 3 == 0:
                        # close out all frame_shifts ending in -1
                        for i in range(len(frame_shifts), 0, -1):
                            if frame_shifts[i - 1][1] < 0:
                                frame_shifts[i - 1][1] = seq[3] + len(seq[0])
                                frame_shifts[i - 1][3] = counter + len(seq[0])
                            else:
                                continue
                        reading_frame = 0
                    else:
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                        reading_frame = (reading_frame + read_frame1 - read_frame2) % 3
            elif seq[2][0][4] == "I":
                # Insertion downstream of start
                coordinates.append(
                    [
                        seq[3],
                        seq[3] + len(seq[0]) * strand - 1,
                        counter,
                        counter + len(seq[0]) - 1,
                        "NA",
                        "NA",
                        seq[2],
                    ]
                )
                if len(seq[0]) % 3 and counter:
                    if not reading_frame:
                        reading_frame = len(seq[0]) % 3
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                    elif not (reading_frame + len(seq[0])) % 3:
                        # close out all frame_shifts ending in -1
                        for i in range(len(frame_shifts), 0, -1):
                            if frame_shifts[i - 1][1] < 0:
                                frame_shifts[i - 1][1] = seq[3] + len(seq[0])
                                frame_shifts[i - 1][3] = counter + len(seq[0])
                            else:
                                continue
                        reading_frame = 0
                    else:
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                        reading_frame = (reading_frame + len(seq[0])) % 3
            # handle a collection of one or more single nucleotide variants
            elif seq[2][0][4] == "V":
                # only document neopeptides corresponding to missense SNVs
                A1 = 3 * ((counter - coding_start) // 3) + coding_start
                B1 = (
                    3 * ((counter + len(seq[0]) - coding_start - 1) // 3)
                    + coding_start
                    + 3
                )
                A2 = 3 * ((ref_counter - ref_start) // 3) + ref_start
                C = 3 * ((seq[3] - coding_start) // 3) + coding_start
                for i in range(0, B1 - A1, 3):
                    A, A_warnings = seq_to_peptide(
                        sequence[(i + A1) : (i + A1 + 3)],
                        mitochondrial=self.mitochondrial,
                        allow_partial_codons=allow_partial_codons,
                    )
                    B, B_warnings = seq_to_peptide(
                        ref_sequence[(i + A2) : (i + A2 + 3)],
                        mitochondrial=self.mitochondrial,
                        allow_partial_codons=allow_partial_codons,
                    )
                    if A != B:
                        # Missense variant
                        if frame_shifts == []:
                            # No upstream frameshifts to complicate paired normal peptides
                            coordinates.append(
                                [
                                    seq[3],
                                    seq[3] + len(seq[0]) * strand - 1,
                                    counter + i * 3,
                                    counter + i * 3 + 2 - 1,
                                    ref_counter + i * 3,
                                    ref_counter + i * 3 + 2 - 1,
                                    seq[2],
                                ]
                            )
                        else:
                            # Upstream frameshifts complicate paired normal peptides
                            coordinates.append(
                                [
                                    seq[3],
                                    seq[3] + len(seq[0]) * strand - 1,
                                    counter + i * 3,
                                    counter + i * 3 + 2 - 1,
                                    "NA",
                                    "NA",
                                    seq[2],
                                ]
                            )
            elif seq[2][0][4] == "E":
                # only document RNA edits corresponding to missense outputs
                A1 = 3 * ((counter - coding_start) // 3) + coding_start
                B1 = (
                    3 * ((counter + len(seq[0]) - coding_start - 1) // 3)
                    + coding_start
                    + 3
                )
                A2 = 3 * ((ref_counter - ref_start) // 3) + ref_start
                C = 3 * ((seq[3] - coding_start) // 3) + coding_start
                for i in range(0, B1 - A1, 3):
                    A, A_warnings = seq_to_peptide(
                        sequence[(i + A1) : (i + A1 + 3)],
                        mitochondrial=self.mitochondrial,
                        allow_partial_codons=allow_partial_codons,
                    )
                    B, B_warnings = seq_to_peptide(
                        ref_sequence[(i + A2) : (i + A2 + 3)],
                        mitochondrial=self.mitochondrial,
                        allow_partial_codons=allow_partial_codons,
                    )
                    if A != B:
                        # Missense variant
                        if frame_shifts == []:
                            # No upstream frameshifts to complicate paired normal peptides
                            coordinates.append(
                                [
                                    seq[3],
                                    seq[3] + len(seq[0]) * strand - 1,
                                    counter + i * 3,
                                    counter + i * 3 + 2 - 1,
                                    ref_counter + i * 3,
                                    ref_counter + i * 3 + 2 - 1,
                                    seq[2],
                                ]
                            )
                        else:
                            # Upstream frameshifts complicate paired normal peptides
                            coordinates.append(
                                [
                                    seq[3],
                                    seq[3] + len(seq[0]) * strand - 1,
                                    counter + i * 3,
                                    counter + i * 3 + 2 - 1,
                                    "NA",
                                    "NA",
                                    seq[2],
                                ]
                            )
        if reading_frame != 0:
            for i in range(len(frame_shifts), 0, -1):
                if frame_shifts[i - 1][1] < 0:
                    frame_shifts[i - 1][1] = annotated_seq[-1][3] + len(
                        annotated_seq[-1][0]
                    )
                    frame_shifts[i - 1][3] = len(sequence)
                else:
                    break
        try:
            protein, editing_positions, ambiguous_positions, protein_warnings = seq_to_peptide(
                    sequence[start_codon[0] : stop_codon[0]],
                    reverse_strand=False,
                    return_positions=True,
                    mitochondrial=self.mitochondrial,
                    allow_partial_codons=allow_partial_codons
                    )
        except TypeError:
            protein, protein_warnings = seq_to_peptide(
                sequence[start_codon[0] :],
                reverse_strand=False,
                mitochondrial=self.mitochondrial,
                allow_partial_codons=allow_partial_codons,
            )
            editing_positions = []
        try:
           protein_ref, ref_editing_positions, ref_ambiguous_positions, ref_warnings = seq_to_peptide(
                ref_sequence[ref_atg[1] : ref_stop[1]],
                reverse_strand=False, return_positions=True,
                mitochondrial=self.mitochondrial,
                allow_partial_codons=allow_partial_codons
                )
        except TypeError:
            ref_editing_positions = []
            try:
                protein_ref, ref_warnings = seq_to_peptide(
                    ref_sequence[ref_atg[1] :],
                    reverse_strand=False,
                    mitochondrial=self.mitochondrial,
                    allow_partial_codons=allow_partial_codons,
                )
            except TypeError:
                protein_ref = ""
        # Check for unknown amino acids
        if "?" in protein or "?" in protein_ref or "X" in protein or "X" in protein_ref:
            unknown_aa = True
        # Turn transcript warnings to tuple of single string
        if not transcript_warnings and not protein_warnings:
            transcript_warnings = ["NA"]
        else:
            transcript_warnings.extend(protein_warnings)
            transcript_warnings = [";".join(transcript_warnings)]
        # Enumerate peptides
        peptide_seqs = collections.defaultdict(list)
        # get amino acid ranges for kmerization
        for size in range(min_size, max_size + 1):
            epitope_coords = []
            peptides_ref = kmerize_peptide(protein_ref, min_size=size, max_size=size, editing_positions=ref_editing_positions)
            for coords in coordinates:
                if coords[4] != "NA" and same_start_frame:
                    # Get coordinates of paired normal peptide
                    epitope_coords.append(
                        [
                            max(0, ((coords[2] - coding_start) // 3) - size + 1),
                            min(len(protein), ((coords[3] - coding_start) // 3) + size),
                            max(0, ((coords[4] - coding_start) // 3) - size + 1),
                            min(
                                len(protein_ref),
                                ((coords[5] - coding_start) // 3) + size,
                            ),
                            coords[6],
                        ]
                    )
                else:
                    # No valid paired normal peptide
                    epitope_coords.append(
                        [
                            max(0, ((coords[2] - coding_start) // 3) - size + 1),
                            min(len(protein), ((coords[3] - coding_start) // 3) + size),
                            "NA",
                            "NA",
                            coords[6],
                        ]
                    )
            for coords in frame_shifts:
                epitope_coords.append(
                    [
                        max(0, ((coords[2] - coding_start) // 3) - size + 1),
                        min(len(protein), ((coords[3] - coding_start) // 3) + size),
                        "NA",
                        "NA",
                        coords[4],
                    ]
                )
            for coords in epitope_coords:
                peptides = kmerize_peptide(
                    protein[coords[0] : coords[1]], min_size=size, max_size=size, editing_positions=editing_positions
                )
                # make normal peptide + ref peptide pair
                if coords[2] != "NA":
                    paired_peptides = kmerize_peptide(
                        protein_ref[coords[2] : coords[3]], min_size=size, max_size=size, editing_positions=ref_editing_positions
                    )
                    if len(paired_peptides) == len(peptides):
                        peptide_pairs = zip(peptides, paired_peptides)
                    else:
                        peptide_pairs = zip(
                            peptides, ["NA" for j in range(0, len(peptides))]
                        )
                    for pair in peptide_pairs:
                        if pair[0] not in peptides_ref:
                            if len(coords[4]) == 2 and type(coords[4][0]) == list:
                                # Dealing with peptide resulting from hybrid interval
                                data_set = coords[4][0][4]
                            else:
                                # Dealing with regular peptide
                                data_set = coords[4]
                            if len(paired_peptides) == len(peptides):
                                if pair[0][1]==True or pair[1][1]==True:
                                    transcript_warnings[0] = ';'.join([transcript_warnings[0], "rna_editing"])
                                if pair[0][2]==True or pair[1][2]==True:
                                    transcript_warnings[0] = ';'.join([transcript_warnings[0],"ambiguous_inosine_codon"])

                            for mutation_data in data_set:
                                if (
                                    unknown_aa
                                    and "?" in pair[0][0]
                                    or "?" in pair[1][0]
                                    or "X" in pair[0][0]
                                    or "X" in pair[1][0]
                                ):
                                    if self.seleno:
                                        mutation_data = (
                                            mutation_data
                                            + (pair[1][0],)
                                            + (
                                                ";".join(
                                                    [
                                                        transcript_warnings[0],
                                                        "unknown_amino_acid",
                                                        "may_contain_selenocysteine",
                                                    ]
                                                ),
                                            )
                                        )
                                    else:
                                        mutation_data = (
                                            mutation_data
                                            + (pair[1][0],)
                                            + (
                                                ";".join(
                                                    [
                                                        transcript_warnings[0],
                                                        "unknown_amino_acid",
                                                    ]
                                                ),
                                            )
                                        )
                                else:
                                    mutation_data = (
                                        mutation_data + (pair[1][0],) + tuple(transcript_warnings)
                                    )
                                peptide_seqs[pair[0][0]].append(mutation_data)
                            peptide_seqs[pair[0][0]] = list(set(peptide_seqs[pair[0][0]]))
                else:
                    peptides = list(set(peptides).difference(peptides_ref))
                    for pep in peptides:
                        if len(coords[4]) == 2 and type(coords[4][0]) == list:
                            # Dealing with peptide resulting from hybrid interval
                            data_set = coords[4][0][4]
                        else:
                            # Dealing with regular peptide
                            data_set = coords[4]

                        if pep[1]==True:
                            transcript_warnings[0] = ';'.join([transcript_warnings[0], "rna_editing"])

                        if pep[2]==True:
                            transcript_warnings[0] = ';'.join([transcript_warnings[0], "ambiguous_inosine_codon"])
                        for mutation_data in data_set:
                            if unknown_aa and "?" in pep[0] or "X" in pep[0]:
                                if self.seleno:
                                    mutation_data = (
                                        mutation_data
                                        + ("NA",)
                                        + (
                                            ";".join(
                                                [
                                                    transcript_warnings[0],
                                                    "unknown_amino_acid",
                                                    "may_contain_selenocysteine",
                                                ]
                                            ),
                                        )
                                    )
                                else:
                                    mutation_data = (
                                        mutation_data
                                        + ("NA",)
                                        + (
                                            ";".join(
                                                [
                                                    transcript_warnings[0],
                                                    "unknown_amino_acid",
                                                ]
                                            ),
                                        )
                                    )
                            else:
                                mutation_data = (
                                    mutation_data + ("NA",) + tuple(transcript_warnings)
                                )
                            peptide_seqs[pep[0]].append(mutation_data)
                        peptide_seqs[pep[0]] = list(set(peptide_seqs[pep[0]]))
        if not return_protein:
            # return list of unique neoepitope sequences
            return peptide_seqs
        else:
            # return list of unique neoepitope sequences plus whole protein
            return peptide_seqs, protein


def gtf_to_cds(gtf_file, dict_dir=None):
    """References cds_dict to get cds bounds for later Bowtie query
    Keys in the dictionary are transcript IDs, while entries are lists of
        relevant CDS/stop codon data
        Data: [chromosome, sequence type, start, stop,
                +/- strand, transcript type]
    Writes cds_dict as a pickled dictionary
    gtf_file: input gtf file to process
    dict_dir: path to directory to store pickled dicts
    Return value: dictionaries
    """
    cds_dict = collections.defaultdict(list)
    cds_lines = collections.defaultdict(list)
    tx_data_dict = collections.defaultdict(list)
    # Parse GTF to obtain CDS/stop codon info
    with xopen(None, gtf_file) as f:
        for line in f:
            try:
                line = line.decode("ascii")
            except AttributeError:
                # it's a string
                pass
            if line[0] != "#":
                tokens = line.strip().split("\t")
                if tokens[2] in ["exon", "start_codon", "stop_codon"]:
                    transcript_id = re.sub(
                        r".*transcript_id \"([A-Z0-9._]+)\"[;].*", r"\1", tokens[8]
                    )
                    transcript_type = re.sub(
                        r".*transcript_type \"([A-Za-z_]+)\"[;].*", r"\1", tokens[8]
                    )
                    if transcript_type in [
                        "protein_coding",
                        "nonsense_mediated_decay",
                        "polymorphic_pseudogene",
                        "IG_V_gene",
                        "TR_V_gene",
                    ]:
                        # Create new dictionary entry for new transcripts
                        cds_dict[transcript_id].append(
                            [
                                tokens[0],
                                tokens[2],
                                int(tokens[3]),
                                int(tokens[4]),
                                tokens[6],
                                transcript_type,
                            ]
                        )
                elif tokens[2] == "CDS":
                    transcript_id = re.sub(
                        r".*transcript_id \"([A-Z0-9._]+)\"[;].*", r"\1", tokens[8]
                    )
                    cds_lines[transcript_id].append(tokens)
                elif tokens[2] == "transcript":
                    transcript_id = re.sub(
                        r".*transcript_id \"([A-Z0-9._]+)\"[;].*", r"\1", tokens[8]
                    )
                    transcript_type = re.sub(
                        r".*transcript_type \"([A-Za-z_]+)\"[;].*", r"\1", tokens[8]
                    )
                    gene_id = re.sub(
                        r".*gene_id \"([A-Z0-9._]+)\"[;].*", r"\1", tokens[8]
                    )
                    gene_name = re.sub(
                        r".*gene_name \"([A-Za-z0-9._-]+)\"[;].*", r"\1", tokens[8]
                    )
                    tag_expr = re.compile('tag "([A-Za-z0-9_]+)"')
                    tags = tag_expr.findall(tokens[8])
                    tsl = None
                    support_level = re.sub(
                        r".*transcript_support_level \"([0-9A-Z]+)\"[;].*",
                        r"\1",
                        tokens[8],
                    )
                    try:
                        tsl = int(support_level)
                    except ValueError:
                        if support_level == "NA":
                            tsl = "NA"
                    if transcript_type in [
                        "protein_coding",
                        "nonsense_mediated_decay",
                        "polymorphic_pseudogene",
                        "IG_V_gene",
                        "TR_V_gene",
                    ]:
                        tx_data_dict[transcript_id] = [
                            transcript_type,
                            gene_id,
                            gene_name,
                            tags,
                            tsl,
                        ]
    # Sort cds_dict coordinates (left -> right) for each transcript
    delete_txs = set()
    for transcript_id, tx_data in cds_dict.items():
        current_cds = cds_lines[transcript_id]
        cds_dict[transcript_id].sort(key=lambda x: x[0])
        seq_types = [x[1] for x in cds_dict[transcript_id]]
        # Assign start codon
        if "start_codon" not in seq_types:
            # Fake a start codon if we have strand info
            try:
                reverse_strand = current_cds[0][6] == "-"
            except IndexError:
                # Remove incompletely annotated transcript
                delete_txs.add(transcript_id)
            else:
                if reverse_strand:
                    # Get transcript intervals
                    intervals = []
                    for seq in [x for x in cds_dict[transcript_id] if x[1] == "exon"]:
                        intervals.extend([seq[2] - 1, seq[3]])
                    intervals.sort()
                    current_cds.sort(key=lambda x: int(x[4]), reverse=True)
                    pos = int(current_cds[0][4]) - int(current_cds[0][7])
                    if bisect.bisect_left(intervals, pos) == bisect.bisect_left(
                        intervals, pos - 2
                    ):
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "start_codon_faux",
                                pos - 2,
                                pos,
                                "-",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    else:
                        index = bisect.bisect_left(intervals, pos)
                        dif = 2 - (pos - intervals[index - 1])
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "start_codon_faux",
                                intervals[index - 2] - dif,
                                intervals[index - 2],
                                "-",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                else:
                    current_cds.sort(key=lambda x: int(x[3]))
                    pos = int(current_cds[0][3]) + int(current_cds[0][7])
                    cds_dict[transcript_id].append(
                        [
                            current_cds[0][0],
                            "start_codon_faux",
                            pos,
                            pos + 2,
                            "+",
                            cds_dict[transcript_id][0][5],
                        ]
                    )
        else:
            # Use annotated start codon
            start_codon_blocks = [
                block for block in cds_dict[transcript_id] if block[1] == "start_codon"
            ]
            if len(start_codon_blocks) > 1:
                min_start = min([int(block[2]) for block in start_codon_blocks])
                for block in start_codon_blocks:
                    if int(block[2]) != min_start:
                        cds_dict[transcript_id].remove(block)
        if "stop_codon" not in seq_types:
            # Fake a start codon if we have strand info
            try:
                reverse_strand = current_cds[0][6] == "-"
            except IndexError:
                # Remove incompletely annotated transcript
                delete_txs.add(transcript_id)
            else:
                # Get transcript intervals
                intervals = []
                for seq in [x for x in cds_dict[transcript_id] if x[1] == "exon"]:
                    intervals.extend([seq[2] - 1, seq[3]])
                intervals.sort()
                # Determine faux stop
                if reverse_strand:
                    # Reverse strand transcript
                    current_cds.sort(key=lambda x: int(x[3]))
                    cds_hang = (
                        (int(current_cds[0][4]) - int(current_cds[0][3]) + 1)
                        - int(current_cds[0][7])
                    ) % 3
                    pos = int(current_cds[0][3])
                    if cds_hang:
                        pos -= 3 - cds_hang
                    index = bisect.bisect_left(intervals, pos + 2)
                    if (index == bisect.bisect_left(intervals, pos)) and (index % 2):
                        # Start and end of stop codon are in the same exon
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "stop_codon_faux",
                                pos,
                                pos + 2,
                                "-",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    elif (
                        (index != bisect.bisect_left(intervals, pos))
                        and (index % 2)
                        and (pos > intervals[0])
                    ):
                        # Stop codon exists partially in an exon and partially in an intron that occurs before the last exon of the transcript
                        dif = (
                            intervals[index - 1] + 1
                        ) - pos  # number of bases we need from previous exon
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "stop_codon_faux",
                                intervals[index - 2] - dif + 1,
                                intervals[index - 2],
                                "-",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    else:
                        # Stop codon overlaps the end of the transcript, so just translate to end
                        continue
                else:
                    # Forward strand transcript
                    current_cds.sort(key=lambda x: int(x[3]), reverse=True)
                    cds_hang = (
                        (int(current_cds[0][4]) - int(current_cds[0][3]) + 1)
                        - int(current_cds[0][7])
                    ) % 3
                    pos = int(current_cds[0][4]) - 2
                    if cds_hang:
                        pos += 3 - cds_hang
                    index = bisect.bisect_left(intervals, pos)
                    if (index == bisect.bisect_left(intervals, pos + 2)) and (
                        index % 2
                    ):
                        # Start and end of stop codon are in the same exon
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "stop_codon_faux",
                                pos,
                                pos + 2,
                                "+",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    elif (
                        (index == bisect.bisect_left(intervals, pos + 2))
                        and (not index % 2)
                        and (index != len(intervals))
                    ):
                        # Stop codon has been shifted into an intron that occurs before the last exon of the transcript
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "stop_codon_faux",
                                intervals[index + 1] + 1,
                                min(intervals[index + 1] + 3, intervals[index + 2]),
                                "+",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    elif (
                        (index != bisect.bisect_left(intervals, pos + 2))
                        and (index % 2)
                        and (index != len(intervals) - 1)
                    ):
                        # Stop codon exists partially in an exon and partially in an intron that occurs before the last exon of the transcript
                        cds_dict[transcript_id].append(
                            [
                                current_cds[0][0],
                                "stop_codon_faux",
                                pos,
                                intervals[index],
                                "+",
                                cds_dict[transcript_id][0][5],
                            ]
                        )
                    else:
                        # Stop codon overlaps the end of the transcript, so just translate to end
                        continue
        else:
            # Assign stop codon if annotated
            stop_codon_blocks = [
                block for block in cds_dict[transcript_id] if block[1] == "stop_codon"
            ]
            if len(stop_codon_blocks) > 1:
                min_stop = min([int(block[2]) for block in stop_codon_blocks])
                for block in stop_codon_blocks:
                    if int(block[2]) != min_stop:
                        cds_dict[transcript_id].remove(block)
    for transcript_id in delete_txs:
        del cds_dict[transcript_id]
    # Write to pickled dictionary
    if dict_dir is not None:
        pickle_dict = os.path.join(dict_dir, "transcript_to_CDS.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(cds_dict, f)
        pickle_dict2 = os.path.join(dict_dir, "transcript_to_gene_info.pickle")
        with open(pickle_dict2, "wb") as f:
            pickle.dump(tx_data_dict, f)
    return cds_dict, tx_data_dict


def cds_to_feature_length(cds_dict, tx_data_dict, dict_dir=None):
    """Creates a dictionary linking gene ID to gene length
    Gene length is median length of all isoforms

    cds_dict: CDS dictionary produced by gtf_to_cds()
    tx_data_dict: transcript data dictionary produced by gtf_to_cds()
    dict_dir: where to write picked dictionary or None if pickling
        shouldn't happen

    Return value: dictionary, keys are gene IDs, values are gene lengths (kilobase)
    """
    gene_to_isoform_lengths = collections.defaultdict(list)
    feature_to_feature_length = {}
    # Iterate through all transcripts to get isoform lengths
    for transcript_id in cds_dict:
        length = 0
        for block in cds_dict[transcript_id]:
            if block[1] == "exon":
                # Calculate exon length and add to isoform length
                block_length = block[3] - block[2] + 1
                length += block_length
        # Store isoform length in dict
        feature_to_feature_length[transcript_id] = length / 1000.0
    # Write to pickled dictionary
    if dict_dir is not None:
        pickle_dict = os.path.join(dict_dir, "feature_to_feature_length.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(feature_to_feature_length, f)
    return feature_to_feature_length


def cds_to_tree(cds_dict, dict_dir=None):
    """Creates searchable tree of chromosome intervals from CDS dictionary
    
    Each chromosome is stored in the dictionary as an interval tree object
    Intervals are added for each CDS, with the associated transcript ID
    Assumes transcript is all on one chromosome - does not work for gene fusions
    Writes the searchable tree as a pickled dictionary

    cds_dict: CDS dictionary produced by gtf_to_cds()
    dict_dir: where to write picked dictionary or None if pickling
        shouldn't happen

    Return value: searchable tree
    """
    searchable_tree = {}
    # Add genomic intervals to the tree for each transcript
    for transcript_id in cds_dict:
        transcript = cds_dict[transcript_id]
        chrom = transcript[0][0]
        # Add new entry for chromosome if not already encountered
        if chrom not in searchable_tree:
            searchable_tree[chrom] = IntervalTree()
        # Add CDS interval to tree with transcript ID
        for cds in transcript:
            start = cds[2]
            stop = cds[3] + 1
            # Interval coordinates are inclusive of start, exclusive of stop
            if stop > start:
                searchable_tree[chrom][start:stop] = transcript_id
            # else:
            # report an error?
    # Write to pickled dictionary
    if dict_dir is not None:
        pickle_dict = os.path.join(dict_dir, "intervals_to_transcript.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(searchable_tree, f)
    return searchable_tree

def transcript_to_rna_edits(rna_edit_file, cds_tree, cds_dict,
                             dict_dir=None):
    """ Creates dictionary mapping transcripts to RNA edit positions

    rna_edit_file: path to REDIportal-formatted file storing RNA A-to-I
        edits. The first line of this file is a header that can be
        anything, and subsequent lines must have, at a minimum, 5 tab-
        separated fields:
            chrom [TAB] 1-based position [TAB] Reference base (can be
            blank) [TAB] Edit base (can be blank) [TAB] Strand (+/-)
    cds_tree: tree output by cds_to_tree() in this script
    cds_dict: cds_dict output by gtf_to_cds() 
    dict_dir: where to write picked dictionary or None if pickling
        shouldn't happen

    Return value: dictionary mapping each transcript ID to a list of
    tuples recording RNA edits in the format
    (chromosome, 0-based position)
    """
    transcript_to_rna_edits = collections.defaultdict(list)
    with xopen(None, rna_edit_file) as f:
        # remove header from rediportal file
        next(f)
        for line in f:
            chrom, pos, _, _, strand = line.split('\t')[:5]
            pos = int(pos)
            # 0.045% of rediportal region labels are alt chromosomes (e.g. 'chr14_GL000009v2_random')
            try:
                cds_tree[chrom]
                for transcript in cds_tree[chrom][pos:pos + 1]:
                    # some transcript keys returning empty list values from cds_dict
                    if not cds_dict[transcript]:
                        pass
                    elif strand == cds_dict[transcript][0][4]:
                        transcript_to_rna_edits[transcript].append((chrom, pos))
            # alt chroms not in cds_tree skipped
            except KeyError:
                pass
    if dict_dir is not None:
        pickle_dict = os.path.join(dict_dir, "transcript_to_rna_edits.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(transcript_to_rna_edits, f)
    return transcript_to_rna_edits

def get_transcripts_from_tree(chrom, start, stop, cds_tree):
    """Uses cds tree to btain transcript IDs from genomic coordinates

    chrom: (String) Specify chrom to use for transcript search.
    start: (Int) Specify start position to use for transcript search.
    stop: (Int) Specify ending position to use for transcript search
    cds_tree: (Dict) dictionary of IntervalTree() objects containing
        transcript IDs as function of exon coords indexed by chr/contig ID.

    Return value: (set) a set of matching unique transcript IDs.
    """
    transcript_ids = set()
    # Interval coordinates are inclusive of start, exclusive of stop
    if chrom not in cds_tree:
        return []
    cds = list(cds_tree[chrom].overlap(start, stop))
    for cd in cds:
        transcript_ids.add(cd.data)
    return list(transcript_ids)


def add_mut_to_haplotype_block(alternative, contig):
    ## if len(ALTERNATIVE) > 1? plural not in args
    if len(alternatives) > 1:
        if i == 0:
            if tokens[1] == "1":
                gen1 = "1"
                gen2 = "0"
            else:
                gen1 = "0"
                gen2 = "1"
        elif i == 1:
            if tokens[1] == "2":
                gen1 = "1"
                gen2 = "0"
            else:
                gen1 = "0"
                gen2 = "1"
    else:
        gen1 = tokens[1]
        gen2 = tokens[2]
    overlapping_transcripts = get_transcripts_from_tree(contig, pos, end, interval_dict)
    if not phasing or gen1 != gen2:
        # For each overlapping transcript, add mutation entry
        # Contains chromosome, position, reference, alternate, allele
        #   A, allele B, genotype line from VCF
        for transcript in overlapping_transcripts:
            block_transcripts[transcript].append(
                [
                    contig,
                    pos,
                    ref,
                    alt,
                    gen1,
                    gen2,
                    tokens[7],
                    mutation_type,
                ]
            )
    else:
        for transcript in overlapping_transcripts:
            homozygous_variants[transcript].append(
                [
                    contig,
                    pos,
                    ref,
                    alt,
                    gen1,
                    gen2,
                    tokens[7],
                    mutation_type,
                ]
            )


def process_haplotypes(hapcut_output, interval_dict, phasing):
    """Stores all haplotypes relevant to different transcripts as a dictionary
    hapcut_output: output from HAPCUT2, adjusted to include unphased
                    mutations as their own haplotypes (performed in
                    software's prep mode)
    interval_dict: dictionary linking genomic intervals to transcripts
    phasing: whether to phase mutations (boolean)
    Return value: dictinoary linking haplotypes to transcripts
    """
    chr_in_intervals = False
    for contig in interval_dict:
        if "chr" in contig:
            chr_in_intervals = True
            continue
    affected_transcripts = collections.defaultdict(list)
    homozygous_variants = collections.defaultdict(list)
    try:
        if hapcut_output == "-":
            input_stream = sys.stdin
        else:
            input_stream = open(hapcut_output)
        block_transcripts = collections.defaultdict(list)
        block_complex_pairs = []
        for line in input_stream:
            if line.startswith("BLOCK"):
                # Skip block header lines
                continue
            elif line[0] == "*":
                # Process all transcripts for the block
                for transcript_id in block_transcripts:
                    block_transcripts[transcript_id].sort(key=itemgetter(1))
                    if phasing:
                        haplotype = []
                        for mut in block_transcripts[transcript_id]:
                            haplotype.append(mut)
                        affected_transcripts[transcript_id].append(haplotype)
                    else:
                        paired_muts = []
                        # First add mutations broken down from complex indels as haplotypes
                        for pair in block_complex_pairs:
                            affected_transcripts[transcript_id].append(pair)
                            paired_muts.extend(pair)
                        # Then add simple mutations as their own haplotypes
                        for mut in block_transcripts[transcript_id]:
                            if mut not in paired_muts:
                                affected_transcripts[transcript_id].append([mut])
                # Reset transcript dictionary
                block_transcripts = collections.defaultdict(list)
                block_complex_pairs = []
            else:
                # Add mutation to transcript dictionary for the block
                tokens = line.strip("\n").split()
                contig = tokens[3]
                if (
                    chr_in_intervals
                    and "chr" not in contig
                    and "".join(["chr", contig]) in interval_dict
                ):
                    contig = "chr" + contig
                if "," in tokens[6]:
                    alternatives = tokens[6].split(",")
                else:
                    alternatives = [tokens[6]]
                for i in range(0, min(len(alternatives), 2)):
                    variants_to_process = []
                    if alternatives[i] == "<DEL>" or alternatives[i] == "*":
                        mutation_type = "D"
                        pos = int(tokens[4])
                        deletion_size = len(tokens[5])
                        ref = tokens[5]
                        alt = deletion_size
                        end = pos + deletion_size
                        variants_to_process.append([pos, ref, alt, end, mutation_type])
                    elif len(tokens[5]) == len(alternatives[i]):
                        mutation_type = "V"
                        pos = int(tokens[4])
                        ref = tokens[5]
                        alt = alternatives[i]
                        mut_size = len(tokens[5])
                        end = pos + mut_size
                        variants_to_process.append([pos, ref, alt, end, mutation_type])
                    elif len(tokens[5]) > len(alternatives[i]):
                        if tokens[5].startswith(alternatives[i]):
                            # Simple deletion
                            mutation_type = "D"
                            deletion_size = len(tokens[5]) - len(alternatives[i])
                            pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                            ref = tokens[5][len(alternatives[i]) :]
                            alt = deletion_size
                            end = pos + deletion_size
                            variants_to_process.append(
                                [pos, ref, alt, end, mutation_type]
                            )
                        else:
                            # Complex indel
                            # Add deletion first
                            mutation_type = "D"
                            pos = int(tokens[4])
                            ref = tokens[5]
                            alt = len(tokens[5])
                            end = pos + alt
                            variants_to_process.append(
                                [pos, ref, alt, end, mutation_type]
                            )
                            # Then add insertion
                            mutation_type = "I"
                            pos = int(tokens[4]) + len(tokens[5]) - 1
                            ref = ""
                            alt = alternatives[i]
                            end = pos + 1
                            variants_to_process.append(
                                [pos, ref, alt, end, mutation_type]
                            )
                    elif len(tokens[5]) < len(alternatives[i]):
                        if alternatives[i].startswith(tokens[5]):
                            # Simple insertion
                            mutation_type = "I"
                            insertion_size = len(alternatives[i]) - len(tokens[5])
                            pos = int(tokens[4]) + len(tokens[5]) - 1
                            ref = ""
                            alt = alternatives[i][len(tokens[5]) :]
                            end = pos + 1
                            variants_to_process.append(
                                [pos, ref, alt, end, mutation_type]
                            )
                        else:
                            # Complex indel - add deletion first
                            mutation_type = "D"
                            pos = int(tokens[4])
                            ref = tokens[5]
                            alt = len(tokens[5])
                            end = pos + alt
                            variants_to_process.append(
                                [pos, ref, alt, end, mutation_type]
                            )
                            # Then add insertion
                            mutation_type = "I"
                            pos = int(tokens[4]) + len(tokens[5]) - 1
                            ref = ""
                            alt = alternatives[i]
                            end = pos + 1
                            variants_to_process.append(
                                (pos, ref, alt, end, mutation_type)
                            )
                    # Store complex variants together if relevant
                    complex_pairs = []
                    for (pos, ref, alt, end, mutation_type) in variants_to_process:
                        if len(alternatives) > 1:
                            if i == 0:
                                if tokens[1] == "1":
                                    gen1 = "1"
                                    gen2 = "0"
                                else:
                                    gen1 = "0"
                                    gen2 = "1"
                            elif i == 1:
                                if tokens[1] == "2":
                                    gen1 = "1"
                                    gen2 = "0"
                                else:
                                    gen1 = "0"
                                    gen2 = "1"
                        else:
                            gen1 = tokens[1]
                            gen2 = tokens[2]
                        overlapping_transcripts = get_transcripts_from_tree(
                            contig, pos, end, interval_dict
                        )
                        if not phasing or gen1 != gen2:
                            # For each overlapping transcript, add mutation entry
                            # Contains chromosome, position, reference, alternate, allele
                            #   A, allele B, genotype line from VCF
                            for transcript in overlapping_transcripts:
                                block_transcripts[transcript].append(
                                    [
                                        contig,
                                        pos,
                                        ref,
                                        alt,
                                        gen1,
                                        gen2,
                                        tokens[7],
                                        mutation_type,
                                    ]
                                )
                                complex_pairs.append(
                                    [
                                        contig,
                                        pos,
                                        ref,
                                        alt,
                                        gen1,
                                        gen2,
                                        tokens[7],
                                        mutation_type,
                                    ]
                                )
                        else:
                            for transcript in overlapping_transcripts:
                                homozygous_variants[transcript].append(
                                    [
                                        contig,
                                        pos,
                                        ref,
                                        alt,
                                        gen1,
                                        gen2,
                                        tokens[7],
                                        mutation_type,
                                    ]
                                )
                    # Store complex pairs if the variant was complex
                    if len(complex_pairs) > 1:
                        complex_pairs.sort(key=itemgetter(1))
                        block_complex_pairs.append(complex_pairs)
    finally:
        if input_stream is not sys.stdin:
            input_stream.close()
    return affected_transcripts, homozygous_variants


def get_haplotype_cliques(haplotype):
    """Finds the maximal cliques of phased variants for a predicted haplotype.

    HapCUT2 and GATK's ReadBackedPhasing may phase together incompatible,
    overlapping variants. This function turns a predicted haplotype into a
    graph, where variants that are predicted to be phased are connected by
    edges only if they are compatible with each other. The maximal cliques
    are then found and returned as their own haplotypes.

    haplotype: predicted haplotype block (as output from process_haplotypes);
        list of lists containing [chromosome, position, reference allele,
        alternate allele, presence on DNA copy 1 (0/1), presence on DNA copy
        2 (0/1), variant information from VCF, variant type ('V', 'D', or 'I')]
        for each variant in the block

    Return value: list of maximal cliques within the haplotype
    """
    graph = nx.Graph()
    for i in range(len(haplotype)):
        # Add node to graph for variant
        graph.add_node(tuple(haplotype[i]))
        # Check mutation class and start position of variant
        if "*" in haplotype[i][6]:
            i_class = "G"
        else:
            i_class = "S"
        i_start = haplotype[i][1]
        # Iterate through other variants to check compatibility
        for j in range(len(haplotype)):
            # Variants are different and predicted to be phased together
            if j != i and (
                (haplotype[i][4] == haplotype[j][4])
                or (haplotype[i][5] == haplotype[j][5])
            ):
                # Check mutation class and start position of second variant
                if "*" in haplotype[j][6]:
                    j_class = "G"
                else:
                    j_class = "S"
                j_start = haplotype[j][1]
                # Check for compatibility of the two variants, add edge to graph if they're compatible
                if (
                    haplotype[i][7] == "I"
                    and haplotype[j][7] == "I"
                    and haplotype[i][1] == haplotype[j][1]
                ):
                    # Insertions that start at same position are incompatible
                    continue
                elif (
                    haplotype[i][7] == "V"
                    and haplotype[j][7] == "V"
                    and i_class == j_class
                ):
                    # Check whether substitutions of same mutation class overlap
                    i_end = i_start + len(haplotype[i][3]) - 1
                    j_end = j_start + len(haplotype[j][3]) - 1
                    if (i_end >= j_start and i_end <= j_end) or (
                        i_start >= j_start and i_start <= j_end
                    ):
                        # Substitutions overlap and are incompatible
                        continue
                    elif (i_start >= j_start and i_end <= j_end) or (
                        j_start >= i_start and j_end <= i_end
                    ):
                        # One substitution completely inside the other
                        continue
                    else:
                        # Substitutions are compatible
                        graph.add_edge(tuple(haplotype[i]), tuple(haplotype[j]))
                elif (
                    haplotype[i][7] == "D"
                    and haplotype[j][7] == "D"
                    and i_class == j_class
                ):
                    # Check whether deletions of same mutation class overlap
                    i_end = i_start + haplotype[i][3] - 1
                    j_end = j_start + haplotype[j][3] - 1
                    if (i_end >= j_start and i_end <= j_end) or (
                        i_start >= j_start and i_start <= j_end
                    ):
                        # Deletions overlap and are incompatible
                        continue
                    elif (i_start >= j_start and i_end <= j_end) or (
                        j_start >= i_start and j_end <= i_end
                    ):
                        # One deletion completely inside the other
                        continue
                    else:
                        # Deletions are compatible
                        graph.add_edge(tuple(haplotype[i]), tuple(haplotype[j]))
                else:
                    # Variants are compatible
                    graph.add_edge(tuple(haplotype[i]), tuple(haplotype[j]))
    # Return maximal cliques for this predicted haplotype
    return list(nx.find_cliques(graph))


def get_peptides_from_transcripts(
    relevant_transcripts,
    homozygous_variants,
    vaf_pos,
    cds_dict,
    info_dict,
    only_novel_upstream,
    only_downstream,
    only_reference,
    reference_index,
    size_list,
    nmd,
    pp,
    igv,
    trv,
    allow_nonstart,
    allow_nonstop,
    allow_partial_codons=False,
    include_germline=2,
    include_somatic=1,
    include_rna_edits=0,
    protein_fasta=False,
    rna_edit_dict=None
):
    """For transcripts that are affected by a mutation, mutations are applied
    and neoepitopes resulting from mutations are called

    relevant_transcripts: dictionary linking haplotypes to transcripts;
        output from process_haplotypes()
    vaf_pos: position of VAF in VCF mutation data from HapCUT2
    cds_dict: dictionary linking transcript IDs, to lists of
        relevant CDS/stop codon data; output from gtf_to_cds()
    info_dict: dictionary linking transcript IDs to transcript information
    only_novel_upstream: whether to start translation from novel upstream
        start codons (boolean)
    only_downstream: whether to start translation from only downstream of
        a disrupted canonical start codon (boolean)
    only_reference: whether to start translation only from the canonical
        start codon for a transcript
    reference_index: BowtieIndexReference object for retrieving
        reference genome sequence
    size_list: list of peptide sizes for neoepitope enumeration
    include_germline: 0 = do not include germline mutations,
            1 = exclude germline mutations from reference comparison,
            2 = include germline mutations in both annotated sequence and
            reference comparison
    include_somatic: 0 = do not include somatic mutations,
            1 = exclude somatic mutations from reference comparison,
            2 = include somatic mutations in both annotated sequence and
            reference comparison
    include_rna_edits: 0 = do not include A to I RNA editing,
            1 = exclude RNA edits from reference comparison,
            2 = include RNA edits in both annotated sequence and reference comparison
    nmd: whether to include nonsense mediated decay transcripts (boolean)
    pp: whether to include polymorphic pseudogene transcripts (boolean)
    igv: whether to include IGV transcripts (boolean)
    trv: whether to include TRV transcripts (boolean)
    allow_nonstart: whether to allow transcripts w/o start codons
    allow_nonstop: wheather to allow transcripts w/o stop codons
    allow_partial_codons: attempt to translate partial codons at ends of transcripts
    protein_fasta: wheather to generate full-length protein sequences
        for fasta file
    rna_edit_dict: dictionary linking chrom,pos to interval
    return value: dictionary linking neoepitopes to their associated
        metadata
    """
    neoepitopes = collections.defaultdict(list)
    fasta_entries = collections.defaultdict(set)
    used_homozygous_variants = set()
    for affected_transcript in relevant_transcripts:
        # Filter out NMD, polymorphic pseudogene, IG V, TR V transcripts if relevant
        if cds_dict[affected_transcript][0][5] == "nonsense_mediated_decay" and not nmd:
            continue
        elif cds_dict[affected_transcript][0][5] == "polymorphic_pseudogene" and not pp:
            continue
        elif cds_dict[affected_transcript][0][5] == "IG_V_gene" and not igv:
            continue
        elif cds_dict[affected_transcript][0][5] == "TR_V_gene" and not trv:
            continue
        # Filter out transcripts missing start or stop codons if relevant
        if (
            "start_codon" not in [x[1] for x in cds_dict[affected_transcript]]
            and not allow_nonstart
        ):
            continue
        if (
            "stop_codon" not in [x[1] for x in cds_dict[affected_transcript]]
            and not allow_nonstop
        ):
            continue
        # Create transcript object
        seleno = False
        if "seleno" in info_dict[affected_transcript][3]:
            seleno = True
        transcript_a = Transcript(
                reference_index,
                [
                    [str(chrom), "", seq_type, str(start), str(end), ".", strand]
                    for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[
                        affected_transcript
                    ]
                ],
                affected_transcript,
                seleno,
                rna_edit_dict
            )
        # Iterate over haplotypes associated with this transcript
        haplotypes = relevant_transcripts[affected_transcript]
        for ht in haplotypes:
            # Check for homozygous variants on affected transcript
            if affected_transcript in homozygous_variants:
                for homozygous in homozygous_variants[affected_transcript]:
                    ht.append(homozygous)
                    used_homozygous_variants.add(tuple(homozygous))
            # Find maximal cliques
            cliques = get_haplotype_cliques(ht)
            for c in cliques:
                # Make edits for each mutation
                for mutation in c:
                    # Determine if mutation is somatic or germline
                    if mutation[6][-1] == "*":
                        mutation_class = "G"
                    else:
                        mutation_class = "S"
                    # Determine VAF if available
                    vaf = None
                    if vaf_pos is not None and mutation_class == "S":
                        vaf_entry = mutation[6].strip("*").split(":")[vaf_pos[0]]
                        if "," in vaf_entry:
                            vaf_entry = [x for x in vaf_entry.split(",") if x != "."]
                            if len(vaf_entry) > 0:
                                vaf = sum(
                                    [float(x.strip("%")) for x in vaf_entry]
                                ) / len(vaf_entry)
                                if vaf_pos[1] == "FREQ":
                                    vaf = vaf / 100.0
                        else:
                            if vaf_entry.strip("%") != ".":
                                vaf = float(vaf_entry.strip("%"))
                                if vaf_pos[1] == "FREQ":
                                    vaf = vaf / 100.0
                    # Determine which copies variant exists on & make edits
                    transcript_a.edit(
                        mutation[3],
                        mutation[1],
                        mutation_type=mutation[7],
                        mutation_class=mutation_class,
                        vaf=vaf,
                    )
                # Extract neoepitopes
                peptides, protein = transcript_a.neopeptides(
                    min_size=size_list[0],
                    max_size=size_list[-1],
                    include_somatic=include_somatic,
                    include_germline=include_germline,
                    include_rna_edits=include_rna_edits,
                    only_novel_upstream=only_novel_upstream,
                    only_downstream=only_downstream,
                    only_reference=only_reference,
                    allow_partial_codons=allow_partial_codons,
                    return_protein=True
                )
                # Store neoepitopes and their metadata
                for pep in peptides:
                    for meta_data in peptides[pep]:
                        adj_meta_data = meta_data + (transcript_a.transcript_id,)
                        if adj_meta_data not in neoepitopes[pep]:
                            neoepitopes[pep].append(adj_meta_data)
                if protein_fasta:
                    if len(peptides) > 0 and protein != "":
                        fasta_entries[affected_transcript].add(protein)
                transcript_a.reset(reference=True)
    for transcript in homozygous_variants:
        # Filter out NMD, polymorphic pseudogene, IG V, TR V transcripts if relevant
        if cds_dict[transcript][0][5] == "nonsense_mediated_decay" and not nmd:
            continue
        elif cds_dict[transcript][0][5] == "polymorphic_pseudogene" and not pp:
            continue
        elif cds_dict[transcript][0][5] == "IG_V_gene" and not igv:
            continue
        elif cds_dict[transcript][0][5] == "TR_V_gene" and not trv:
            continue
        seleno = False
        if "seleno" in info_dict[transcript][3]:
            seleno = True
        transcript_a = Transcript(
            reference_index,
            [
                [str(chrom), "", seq_type, str(start), str(end), ".", strand]
                for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[
                    transcript
                ]
            ],
            transcript,
            seleno,
            rna_edit_dict
        )
        for mutation in homozygous_variants[transcript]:
            if tuple(mutation) not in used_homozygous_variants:
                # Determine if mutation is somatic or germline
                if mutation[6][-1] == "*":
                    mutation_class = "G"
                else:
                    mutation_class = "S"
                # Determine VAF if available
                vaf = None
                if vaf_pos is not None and mutation_class == "S":
                    vaf_entry = mutation[6].strip("*").split(":")[vaf_pos[0]]
                    if "," in vaf_entry:
                        vaf_entry = [x for x in vaf_entry.split(",") if x != "."]
                        if len(vaf_entry) > 0:
                            vaf = sum([float(x.strip("%")) for x in vaf_entry]) / len(
                                vaf_entry
                            )
                            if vaf_pos[1] == "FREQ":
                                vaf = vaf / 100.0
                    else:
                        if vaf_entry.strip("%") != ".":
                            vaf = float(vaf_entry.strip("%"))
                            if vaf_pos[1] == "FREQ":
                                vaf = vaf / 100.0
                # Make edits
                transcript_a.edit(
                    mutation[3],
                    mutation[1],
                    mutation_type=mutation[7],
                    mutation_class=mutation_class,
                    vaf=vaf,
                )
                # Extract neoepitopes
                peptides, protein = transcript_a.neopeptides(
                    min_size=size_list[0],
                    max_size=size_list[-1],
                    include_somatic=include_somatic,
                    include_germline=include_germline,
                    include_rna_edits=include_rna_edits,
                    only_novel_upstream=only_novel_upstream,
                    only_downstream=only_downstream,
                    only_reference=only_reference,
                    allow_partial_codons=allow_partial_codons,
                    return_protein=True
                )
                # Store neoepitopes and their metadata
                for pep in peptides:
                    for meta_data in peptides[pep]:
                        adj_meta_data = meta_data + (transcript_a.transcript_id,)
                        if adj_meta_data not in neoepitopes[pep]:
                            neoepitopes[pep].append(adj_meta_data)
                if protein_fasta:
                    if len(peptides) > 0 and protein != "":
                        fasta_entries[transcript].add(protein)
            transcript_a.reset(reference=True)
    return neoepitopes, fasta_entries
