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
import pickle
from intervaltree import Interval, IntervalTree
from operator import itemgetter
import sys
import warnings
import contextlib

from sys import version_info

if version_info[0] < 3:
    from string import maketrans
    revcomp_translation_table = maketrans("ATCG", "TAGC")
else:
    revcomp_translation_table = str.maketrans("ATCG", "TAGC")


@contextlib.contextmanager
def xopen(gzipped, *args):
    """ Passes args on to the appropriate opener, gzip or regular.
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
            with open(args[0], 'rb') as binary_input_stream:
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
                #old_read_eof = gzip.GzipFile._read_eof
                #gzip.GzipFile._read_eof = lambda *args, **kwargs: None
                fh = gzip.open(*args)
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


def custom_bisect_left(a, x, lo=0, hi=None, getter=0):
    """ Same as bisect.bisect_left, but compares only index "getter"
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


def kmerize_peptide(peptide, min_size=8, max_size=11):
    """ Obtains subsequences of a peptide.
        normal_peptide: normal peptide seq
        min_size: minimum subsequence size
        max_size: maximum subsequence size
        Return value: list of all possible subsequences of size between
            min_size and max_size
    """
    peptide_size = len(peptide)
    return [
        item
        for sublist in [
            [peptide[i : i + size] for i in range(peptide_size - size + 1)]
            for size in range(min_size, max_size + 1)
        ]
        for item in sublist
        if "X" not in item
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


def seq_to_peptide(seq, reverse_strand=False, require_ATG=False):
    """ Translates nucleotide sequence into peptide sequence.
        All codons including and after stop codon are recorded as X's.
        seq: nucleotide sequence
        reverse_strand: True iff strand is -
        require_ATG: True iff search for start codon (ATG)
        Return value: peptide string
    """
    if reverse_strand:
        seq = seq[::-1].translate(_complement_table)
    if require_ATG:
        start = seq.find("ATG")
        if start >= 0:
            seq = seq[start:]
        else:
            return ""
    seq_size = len(seq)
    peptide = []
    for i in range(0, seq_size - seq_size % 3, 3):
        if 'N' not in seq[i : i + 3]:
            codon = _codon_table[seq[i : i + 3]]
        elif seq[i : i + 3].count('N') == 1 and seq[i+2] == 'N':
            # Only 1 N in the wobble position
            codon_options = set(
                [_codon_table[''.join([seq[i : i + 2], x])] 
                                        for x in ['A', 'C', 'G', 'T']]
                )
            if len(codon_options) == 1:
                codon = list(codon_options)[0]
            else:
                codon = '?'
        else:
            codon = '?'
        peptide.append(codon)
        if codon == "X":
            break
    # for j in range(i + 3, seq_size - seq_size % 3, 3):
    # peptide.append('X')
    return "".join(peptide)


class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype. """

    # Should we handle somatic deletions that overlap germline mutations?
    # I.E., should we break up a somatic deletion into two separate mutations
    # that surround the germline mutation? Or do we only call the somatic?

    def __init__(self, bowtie_reference_index, cds, transcript_id):
        """ Initializes Transcript object.
            This class assumes edits added to a transcript are properly
            phased, consistent, and nonredundant. Most conspicuously, there
            shouldn't be SNVs or insertions among deleted bases.
            bowtie_reference_index: BowtieIndexReference object for retrieving
                reference genome sequence
            cds: list of all CDS lines for exactly one transcript from GTF;
                a line can be a list pre-split by '\t' or not yet split
            transcript_id: transcript ID
        """
        assert len(cds) > 0
        self.bowtie_reference_index = bowtie_reference_index
        self.transcript_id = transcript_id
        self.intervals = []
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
                if line[2].endswith("_faux"):
                    self._find_new_start = False
                else:
                    # 1-based public, 0-based private
                    self._find_new_start = True
                self.start_codon = int(line[3])
                self._start_codon = self.start_codon - 1
            elif line[2] == "stop_codon":
                self.stop_codon = int(line[3])
                self._stop_codon = self.stop_codon - 1
            else:
                raise NotImplementedError("GTF sequence type not currently supported")
            last_chrom, last_strand = line[0], line[6]
        # Store edits to coding sequence only
        self.edits = collections.defaultdict(list)
        self.deletion_intervals = []
        self.chrom = last_chrom
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
        if self._start_codon:
            self.start_codon_index = bisect.bisect_left(
                self.intervals, self._start_codon
            )
        else:
            self.start_codon_index = None
        if self.stop_codon:
            self.stop_codon_index = bisect.bisect_left(self.intervals, self._stop_codon)
        else:
            self.stop_codon_index = None

    def reset(self, reference=False):
        """ Resets to last save point or reference (i.e., removes all edits).
            reference: if False, tries to reset to last save point, and if that
                doesn't exist, resets to reference. If True, resets to
                reference.
            No return value.
        """
        if reference:
            self.edits = collections.defaultdict(list)
            self.deletion_intervals = []
        else:
            self.edits = copy.copy(self.last_edits)
            self.deletion_intervals = copy.copy(self.last_deletion_intervals)

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
            mutation_type: V for SNV, I for insertion, D for deletion
            mutation_class: S for somatic, G for germline
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
                if seq == ref_deletion:
                    del_interval = (
                        pos - 2,
                        pos + deletion_size - 2,
                        mutation_class,
                        (self.chrom, pos, seq, "", mutation_type, vaf),
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
                    ),
                )
            for interval in self.deletion_intervals:
                if del_interval[2] == interval[2]:
                    if (
                        del_interval[0] >= interval[0]
                        and del_interval[0] <= interval[1]
                    ) or (
                        del_interval[1] >= interval[0]
                        and del_interval[1] <= interval[1]
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
            self.edits[pos - 1].append(
                (
                    seq,
                    mutation_type,
                    mutation_class,
                    (self.chrom, pos, "", seq, mutation_type, vaf),
                )
            )
        elif mutation_type == "V":
            reference_seq = self.bowtie_reference_index.get_stretch(
                self.chrom, pos - 1, len(seq)
            )
            other_snvs = [edit for edit in self.edits[pos - 1]]
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
        else:
            raise NotImplementedError("Mutation type not yet implemented")

    def expressed_edits(
        self, start=None, end=None, genome=True, include_somatic=1, include_germline=2
    ):
        """ Gets expressed set of edits and transcript intervals.
            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified
            include_somatic: whether to include somatic mutations (boolean)
            include_germline: whether to include germline mutations (boolean)
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
                        relevant_deletion_intervals.extend([
                                    (intervals[i], "R", tuple())
                                    for i in range(start_index, end_index)
                                ])
                    else:
                        if (
                            start_index % 2
                            #or deletion_intervals[i][0] == intervals[start_index]
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
            + relevant_deletion_intervals
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
                    and edit[2] == "S"
                    or include_germline
                    and edit[2] == "G"
                ):
                    if edit[1] == "V":
                        if start_index % 2 and edit[3][1] != edit[0]:
                            # Add edit if and only if it lies within bounds
                            edits[pos].append(edit)
                    elif edit[1] == "I":
                        try:
                            if start_index % 2 or pos == intervals[start_index][0]:
                                # An insertion is valid before or after a block
                                edits[pos].append(edit)
                        except IndexError:
                            continue
            # If there is more than 1 SNV at the same position, one must be
            # germline and the other somatic, as only 1 mutation per mutation
            # class is allowed at the same position. Favor somatic mutation.
            if pos in edits:
                snvs = [v for v in edits[pos] if v[1] == "V"]
                if len(snvs) > 1:
                    germ = [v for v in snvs if v[2] == "G"][0]
                    edits[pos].remove(germ)
        return (edits, adjusted_intervals)

    def save(self):
        """ Creates save point for edits.

            No return value.
        """
        self.last_edits = copy.copy(self.edits)
        self.last_deletion_intervals = copy.copy(self.deletion_intervals)

    def reading_frame(self, pos):
        """ Retrieves reading frame (0, 1, or 2) at given coordinate.

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
            or not self.stop_codon_index
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
                if pos > self._start_codon or pos < self._stop_codon:
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
                if pos < self._start_codon or pos > self._stop_codon:
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
        """ Appends mutation to seq_list, merging successive mutations.
            seq_list: list of tuples (sequence, type) where type is one
                of R, G, or S (for respectively reference, germline edit, or
                somatic edit). Empty sequence means there was a deletion.
            seq: seq to add
            mutation_class: S for somatic, G for germline, R for reference
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
                        [mutation_info[i] for i in range(0, len(mutation_info))],
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
        adj_alt_seq = copy.copy(alt[2])
        adj_alt_allele = copy.copy(alt[3])
        adj_alt_mut_info = copy.copy(alt[4])
        adj_ref_seq = copy.copy(ref[2])
        adj_ref_allele = copy.copy(ref[3])
        adj_ref_mut_info = copy.copy(ref[4])
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
        alt[2] = copy.copy(adj_alt_seq)
        alt[3] = copy.copy(adj_alt_allele)
        alt[4] = copy.copy(adj_alt_mut_info)
        ref[2] = copy.copy(adj_ref_seq)
        ref[3] = copy.copy(adj_ref_allele)
        ref[4] = copy.copy(adj_ref_mut_info)
        return ("", "H", [alt, ref], alt[1])

    def annotated_seq(
        self, start=None, end=None, genome=True, include_somatic=1, include_germline=2
    ):
        """ Retrieves transcript sequence between start and end coordinates.
            Includes info on whether edits are somatic or germline and whether
            sequence is reference sequence.
            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified
            include_somatic: whether to include somatic mutations (boolean)
            include_germline: whether to include germline mutations (boolean)
            Return value: list of tuples (sequence, mutation class,
                mutation information, position),
                where sequence is a segment of sequence of the (possibly)
                mutated transcript, mutation class is one of {'G', 'S', 'H', 'R'},
                where 'G' denotes germline, 'S' denotes somatic, 'H' denotes hybrid
                of somatic and germline, and 'R' denotes reference sequence,
                mutation information is a list of tuple(s) (chromosome,
                1-based position of {first base of deletion, base before
                insertion, SNV}, reference sequence, variant sequence,
                {'D', 'I', 'V'}, VAF) , and position is the 1-based position
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
            )
            """Check for insertions at beginnings of intervals, and if they're
            present, shift them to ends of previous intervals so they're
            actually added."""
            new_edits = copy.copy(edits)
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
                    last_index, last_pos = 0, intervals[i - 1][0] + 1
                    for pos_to_add in pos_group:
                        fill = pos_to_add - last_pos
                        if intervals[i - 1][1] != "R":
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
                        if self.rev_strand:
                            self._seq_append(
                                final_seq,
                                seqs[(i - 1) // 2][0][last_index : last_index + fill],
                                "R",
                                tuple(),
                                seqs[(i - 1) // 2][1][0] + last_index + fill - 1,
                                merge=False,
                            )
                        else:
                            self._seq_append(
                                final_seq,
                                seqs[(i - 1) // 2][0][last_index : last_index + fill],
                                "R",
                                tuple(),
                                seqs[(i - 1) // 2][1][0] + last_index,
                                merge=False,
                            )
                        # If no edits, snv is reference and no insertion
                        try:
                            snv = (
                                seqs[(i - 1) // 2][0][last_index + fill],
                                "R",
                                tuple(),
                                seqs[(i - 1) // 2][1][0] + fill,
                            )
                        except IndexError:
                            """Should happen only for insertions at beginning
                            of sequence."""
                            assert (i - 1) // 2 == 0 and not seqs[0][0]
                            snv = ("", "R", tuple(), seqs[(i - 1) // 2][1][0] + fill)
                        insertion = ("", "R", tuple(), seqs[(i - 1) // 2][1][0] + fill)
                        for edit in new_edits[pos_to_add]:
                            if edit[1] == "V":
                                snv = (edit[0], edit[2], edit[3], edit[3][1])
                            else:
                                assert edit[1] == "I"
                                insertion = (edit[0], edit[2], edit[3], edit[3][1])
                        self._seq_append(final_seq, *snv, merge=False)
                        self._seq_append(final_seq, *insertion, merge=False)
                        last_index += fill + 1
                        last_pos += fill + 1
                    if intervals[i - 1][1] != "R":
                        if isinstance(intervals[i - 1][2], list):
                            genomic_position = min([x[1] for x in intervals[i - 1][2]])
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
                    i += 2
                    try:
                        while pos > intervals[i][0]:
                            if intervals[i - 1][1] != "R":
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
                            i += 2
                    except IndexError:
                        if i > len(intervals) - 1:
                            # Done enumerating sequence
                            break
                    pos_group = [pos]
                else:
                    pos_group.append(pos)
            if self.rev_strand:
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
            return adj_seq
        raise NotImplementedError(
            "Retrieving sequence with transcript coordinates not "
            "yet fully supported."
        )

    def _atg_choice(
        self,
        atgs,
        only_novel_upstream=False,
        only_downstream=True,
        only_reference=False,
    ):
        """ Chooses best start codon from a list of start codons.
            atgs: list of lists, each representing a start codon, of the form
                [pos in sequence or -1 if absent, pos in reference seq or -1 if
                 absent, mutation information,
                 True iff downstream of or at start codon,
                 True iff ATG de novo,
                 True iff ATG absent in annotated seq]
            only_novel_upstream: True if and only if only novel start codons
                upstream of original start codon following reference ATGs
                in 5'UTR are to be considered; otherwise, all upstream ATGs
                are considered.
            only_downstream: True if and only if only downstream alternative
                start codons are to be considered; overrides
                only_novel_upstream.
            only_reference: True if and only if only the original start codon
                is allowed; overrides only_downstream and only_novel_upstream
            Return value: tuple (reading_frame, chosen atg from atgs,
                            reference atg, warnings about start codon choice)
        """
        encountered_true_start = False
        atg_priority_list = []
        start_disrupting_muts = []
        start_warnings = []
        for atg in atgs:
            if not encountered_true_start and atg[3] and not atg[4]:
                if not atg[5] and atg_priority_list == []:
                    """If the original start codon is maintained in the edited
                        transcript, immediately return it."""
                    return 0, atg, atg, []
                else:
                    """If the original start codon is missing in the edited
                        transcript, maintain it to keep track of reading frame
                        changes for new start"""
                    coding_ref_start = atg
                    if atg[5]:
                        for seq in atg[2]:
                            for x in seq[2]:
                                start_disrupting_muts.append(x)
                        start_warnings.append("reference_start_codon_disrupted")
            if atg[3] and not atg[4]:
                encountered_true_start = True
            if not only_reference:
                """If not true start codon and accepting non-reference starts,
                    assess the start codon option"""
                if only_downstream:
                    if atg[3] and not atg[5]:
                        atg_priority_list.append(atg)
                elif only_novel_upstream:
                    if (atg[3] or atg[4]) and not atg[5]:
                        atg_priority_list.append(atg)
                elif not only_novel_upstream and not only_downstream:
                    atg_priority_list.append(atg)
        if atg_priority_list == []:
            # No valid ATGs
            return None, None, None, start_warnings
        else:
            # The start codon is the first of the valid start codons
            start_codon = atg_priority_list[0]
            if type(start_codon[2]) == list:
                for i in range(0, len(start_codon[2])):
                    if start_codon[2][i][2] == [()]:
                        start_codon[2][i] = (
                            start_codon[2][i][0],
                            start_codon[2][i][1],
                            start_disrupting_muts,
                            start_codon[2][i][3],
                        )
            else:
                if start_codon[2][2] == [()]:
                    start_codon[2] = [
                        (
                            start_codon[2][0],
                            start_codon[2][1],
                            start_disrupting_muts,
                            start_codon[2][3],
                        )
                    ]
                else:
                    start_codon[2] = [start_codon[2]]
            distance = start_codon[0] - coding_ref_start[0]
            if start_codon[3]:
                direction = "downstream"
            else:
                direction = "upstream"
            if start_codon[4]:
                novelty = "novel"
            else:
                novelty = "preexisting"
            warnings.warn(
                " ".join(
                    [
                        "Using",
                        novelty,
                        direction,
                        "start codon",
                        str(distance),
                        "bp from reference start codon for",
                        self.transcript_id,
                    ]
                )
            )
            start_warnings.append(
                "_".join(
                    [
                        novelty,
                        direction,
                        "start",
                        "codon",
                        str(distance),
                        "bp",
                        "from",
                        "reference",
                        "start",
                        "codon",
                    ]
                )
            )
        return (
            (start_codon[0] - coding_ref_start[0]) % 3,
            start_codon,
            coding_ref_start,
            start_warnings,
        )

    def neopeptides(
        self,
        min_size=8,
        max_size=11,
        include_somatic=1,
        include_germline=2,
        only_novel_upstream=False,
        only_downstream=True,
        only_reference=False,
        return_protein=False,
    ):
        """ Retrieves dict of predicted peptide fragments from transcript that
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
            Return value: dict of peptides of desired length(s) [KEYS] with
                values equivalent to a list of causal variants [VALUES].
        """
        # if no edits to process, then skip all next steps and return {}
        if include_somatic == include_germline and include_somatic != 1:
            if not return_protein:
                return {}
            else:
                return {}, ""
        if not self.edits and not self.deletion_intervals:
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
        annotated_seq = self.annotated_seq(
            include_somatic=include_somatic, include_germline=include_germline
        )
        sequence, ref_sequence = "", ""  # hold flattened nucleotide sequence
        start = self.start_codon  # redundant var; can change
        stop = self.stop_codon  # redundant var; can change
        unknown_aa = False
        if start is None:
            if not return_protein:
                return {}
            else:
                return {}, ""
        # +1 is + strand, -1 is - strand
        strand = 1 - self.rev_strand * 2
        # hold list of ATGs (from 5' UTR, start, and one downstream of start)
        # ATGs structure is: [pos in sequence (-1 if absent), pos in ref seq
        # (-1 if absent), mutation information, is downstream of start codon?,
        # is ATG new in annotated seq?, is ATG missing in annotated seq?]
        ATGs, TAA_TGA_TAG = [], []
        ATG_counter1, ATG_counter2 = 0, 0
        ATG_limit = 2
        coding_start, ref_start, coding_stop, ref_stop = -1, -1, -1, -1
        counter, ref_counter = 0, 0  # hold edited transcript level coordinates
        seq_previous = []
        new_ATG_upstream = False
        transcript_warnings = []
        annotated_seq.append([])
        for seq in annotated_seq:
            if self._find_new_start:
                # build pairwise list of 'ATG's from annotated_seq and reference
                ATG1 = sequence.find("ATG", max(0, ATG_counter1 - 2))
                ATG2 = ref_sequence.find("ATG", max(0, ATG_counter2 - 2))
                ATG_temp1 = ATG_counter1
                ATG_temp2 = ATG_counter2
                while (ATG1 >= 0 or ATG2 >= 0) and ATG_limit > 0:
                    if seq_previous[-1][3] * strand > start * strand:
                        ATG_limit -= 1
                    if ATG1 - ATG_temp1 == ATG2 - ATG_temp2:
                        # Reference and new sequence contain start in same place
                        if (len(sequence) - len(seq_previous[-1][0]) - 1) < ATG1:
                            relevant_seq = [seq_previous[-1]]
                        else:
                            found_all_variants = False
                            relevant_seq = []
                            i = -1
                            while not found_all_variants and i > -1 * len(seq_previous):
                                prev_seq = "".join([x[0] for x in seq_previous[i:]])
                                if (
                                    len(sequence) - len(prev_seq)
                                ) >= ATG1 and seq_previous[i][2] != [()]:
                                    relevant_seq.append(seq_previous[i])
                                if (
                                    len(sequence) - len(seq_previous[-1][0]) - 1
                                ) < ATG1:
                                    found_all_variants = True
                                i -= 1
                            if relevant_seq == []:
                                relevant_seq = [seq_previous[-1]]
                        ATGs.append(
                            [
                                ATG1,
                                ATG2,
                                relevant_seq,
                                ATG2 >= ref_start and ref_start >= 0,
                                False,
                                False,
                            ]
                        )
                        ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                        ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                        ATG1 = sequence.find("ATG", ATG_counter1)
                        ATG2 = ref_sequence.find("ATG", ATG_counter2)
                    elif ATG1 >= 0 and ATG2 < 0:
                        # New sequence contains start codon while reference doesn't
                        if (len(sequence) - len(seq_previous[-1][0]) - 1) < ATG1:
                            relevant_seq = [seq_previous[-1]]
                        else:
                            found_all_variants = False
                            if seq_previous[-1][2] != [()]:
                                relevant_seq = [seq_previous[-1]]
                            else:
                                relevant_seq = []
                            i = -2
                            while not found_all_variants and i > -1 * len(seq_previous):
                                prev_seq = "".join([x[0] for x in seq_previous[i:]])
                                if (
                                    len(sequence) - len(prev_seq)
                                ) >= ATG1 and seq_previous[i][2] != [()]:
                                    relevant_seq.append(seq_previous[i])
                                if (
                                    len(sequence) - len(seq_previous[-1][0]) - 1
                                ) < ATG1:
                                    found_all_variants = True
                                i -= 1
                        ATGs.append(
                            [
                                ATG1,
                                ATG1 - ATG_temp1 + ATG_temp2,
                                relevant_seq,
                                ATG1 >= coding_start and coding_start >= 0,
                                True,
                                False,
                            ]
                        )
                        ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                        ATG1 = sequence.find("ATG", ATG_counter1)
                    elif ATG1 < 0 and ATG2 >= 0:
                        # Reference contains start codon but new sequence doesn't
                        if (len(ref_sequence) - len(seq_previous[-1][0]) - 1) < ATG2:
                            relevant_seq = [seq_previous[-1]]
                        else:
                            found_all_variants = False
                            if seq_previous[-1][2] != [()]:
                                relevant_seq = [seq_previous[-1]]
                            else:
                                relevant_seq = []
                            i = -2
                            while not found_all_variants and i > -1 * len(seq_previous):
                                prev_seq = "".join([x[0] for x in seq_previous[i:]])
                                if (
                                    len(ref_sequence) - len(prev_seq)
                                ) >= ATG2 and seq_previous[i][2] != [()]:
                                    relevant_seq.append(seq_previous[i])
                                if (
                                    len(ref_sequence) - len(seq_previous[-1][0]) - 1
                                ) < ATG2:
                                    found_all_variants = True
                                i -= 1
                        ATGs.append(
                            [
                                ATG2 - ATG_temp2 + ATG_temp1,
                                ATG2,
                                relevant_seq,
                                ATG2 >= ref_start and ref_start >= 0,
                                False,
                                True,
                            ]
                        )
                        ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                        ATG2 = ref_sequence.find("ATG", ATG_counter2)
                    elif ATG1 - ATG_temp1 < ATG2 - ATG_temp2:
                        # Start codon happens in new sequence before ref sequence
                        if (len(sequence) - len(seq_previous[-1][0]) - 1) < ATG1:
                            relevant_seq = [seq_previous[-1]]
                        else:
                            found_all_variants = False
                            if seq_previous[-1][2] != [()]:
                                relevant_seq = [seq_previous[-1]]
                            else:
                                relevant_seq = []
                            i = -2
                            while not found_all_variants and i > -1 * len(seq_previous):
                                prev_seq = "".join([x[0] for x in seq_previous[i:]])
                                if (
                                    len(sequence) - len(prev_seq)
                                ) >= ATG1 and seq_previous[i][2] != [()]:
                                    relevant_seq.append(seq_previous[i])
                                if (
                                    len(sequence) - len(seq_previous[-1][0]) - 1
                                ) < ATG1:
                                    found_all_variants = True
                                i -= 1
                        ATGs.append(
                            [
                                ATG1,
                                ATG1 - ATG_temp1 + ATG_temp2,
                                relevant_seq,
                                ATG1 >= coding_start and coding_start >= 0,
                                True,
                                False,
                            ]
                        )
                        ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                        ATG1 = sequence.find("ATG", ATG_counter1)
                    else:
                        # Start codon happens in ref sequence before new sequence
                        if (len(ref_sequence) - len(seq_previous[-1][0]) - 1) < ATG2:
                            relevant_seq = [seq_previous[-1]]
                        else:
                            found_all_variants = False
                            if seq_previous[-1][2] != [()]:
                                relevant_seq = [seq_previous[-1]]
                            else:
                                relevant_seq = []
                            i = -2
                            while not found_all_variants and i > -1 * len(seq_previous):
                                prev_seq = "".join([x[0] for x in seq_previous[i:]])
                                if (
                                    len(ref_sequence) - len(prev_seq)
                                ) >= ATG2 and seq_previous[i][2] != [()]:
                                    relevant_seq.append(seq_previous[i])
                                if (
                                    len(ref_sequence) - len(seq_previous[-1][0]) - 1
                                ) < ATG2:
                                    found_all_variants = True
                                i -= 1
                        ATGs.append(
                            [
                                ATG2 - ATG_temp2 + ATG_temp1,
                                ATG2,
                                relevant_seq,
                                ATG2 >= ref_start and ref_start >= 0,
                                False,
                                True,
                            ]
                        )
                        ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                        ATG2 = ref_sequence.find("ATG", ATG_counter2)
                ATG_counter1 = len(sequence)
                ATG_counter2 = len(ref_sequence)
            if seq == []:
                break
            seq_previous.append(seq)
            # find transcript-relative coordinates of start codons
            # flatten strings from annotated and reference seqs
            if seq[1] == "R":
                if ref_start < 0 and seq[3] * strand + len(seq[0]) > start * strand:
                    coding_start = (
                        counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                    ref_start = (
                        ref_counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                if stop is not None:
                    if ref_stop < 0 and seq[3] * strand + len(seq[0]) > stop * strand:
                        coding_stop = (
                            counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        ref_stop = (
                            ref_counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                sequence += seq[0]
                ref_sequence += seq[0]
                counter += len(seq[0])
                ref_counter += len(seq[0])
                continue
            elif seq[1] == "H":
                if (
                    coding_start < 0
                    and seq[2][0][1] * strand + len(seq[2][0][2]) > start * strand
                ):
                    coding_start = (
                        counter + (start - seq[2][0][1] + 2 * self.rev_strand) * strand
                    )
                if (
                    ref_start < 0
                    and seq[2][1][1] * strand + len(seq[2][1][2]) > start * strand
                ):
                    ref_start = (
                        ref_counter
                        + (start - seq[2][1][1] + 2 * self.rev_strand) * strand
                    )
                if stop is not None:
                    if (
                        ref_stop < 0
                        and seq[2][1][1] * strand + len(seq[2][1][2]) > stop * strand
                    ):
                        coding_stop = (
                            counter
                            + (stop - seq[2][0][1] + 2 * self.rev_strand) * strand
                        )
                        ref_stop = (
                            ref_counter
                            + (stop - seq[2][1][1] + 2 * self.rev_strand) * strand
                        )
                        TAA_TGA_TAG = [seq[0], seq[1], seq[2][0][4], seq[3]]
                sequence += seq[2][0][3]
                counter += len(seq[2][0][3])
                if self.rev_strand:
                    ref_sequence += seq[2][1][3][::-1].translate(
                        revcomp_translation_table
                    )
                else:
                    ref_sequence += seq[2][1][3]
                ref_counter += len(seq[2][1][3])
                continue
            elif seq[2][0][4] == "D":
                if (
                    ref_start < 0
                    and seq[3] * strand + len(seq[2][0][2]) > start * strand
                ):
                    coding_start = (
                        counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                    ref_start = (
                        ref_counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                if stop is not None:
                    if ref_stop < 0 and seq[3] * strand + len(seq[0]) > stop * strand:
                        coding_stop = (
                            counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        ref_stop = (
                            ref_counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        TAA_TGA_TAG = seq
                if not (
                    (seq[1] == "G" and include_germline == 2)
                    or (seq[1] == "S" and include_somatic == 2)
                ):
                    #                    ref_sequence += seq[0]
                    #                    ref_counter += len(seq[0])
                    #                else:
                    if self.rev_strand:
                        for i in seq[2]:
                            ref_sequence += i[2][::-1].translate(
                                revcomp_translation_table
                            )
                            ref_counter += len(i[2])
                    else:
                        for i in seq[2]:
                            ref_sequence += i[2]
                            ref_counter += len(i[2])
                continue
            elif seq[2][0][4] == "I":
                if ref_start < 0 and seq[3] * strand + len(seq[0]) > start * strand:
                    coding_start = (
                        counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                    ref_start = (
                        ref_counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                if stop is not None:
                    if ref_stop < 0 and seq[3] * strand + len(seq[0]) > stop * strand:
                        coding_stop = (
                            counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        ref_stop = (
                            ref_counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        TAA_TGA_TAG = seq
                sequence += seq[0]
                counter += len(seq[0])
                if (seq[1] == "G" and include_germline == 2) or (
                    seq[1] == "S" and include_somatic == 2
                ):
                    ref_sequence += seq[0]
                    ref_counter += len(seq[0])
                continue
            elif seq[2][0][4] == "V":
                if ref_start < 0 and seq[3] * strand + len(seq[0]) > start * strand:
                    coding_start = (
                        counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                    ref_start = (
                        ref_counter + (start - seq[3] + 2 * self.rev_strand) * strand
                    )
                if stop is not None:
                    if ref_stop < 0 and seq[3] * strand + len(seq[0]) > stop * strand:
                        coding_stop = (
                            counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        ref_stop = (
                            ref_counter + (stop - seq[3] + 2 * self.rev_strand) * strand
                        )
                        TAA_TGA_TAG = seq
                sequence += seq[0]
                counter += len(seq[0])
                if (seq[1] == "G" and include_germline == 2) or (
                    seq[1] == "S" and include_somatic == 2
                ):
                    ref_sequence += seq[0]
                else:
                    if self.rev_strand:
                        for i in seq[2]:
                            ref_sequence += i[2][::-1].translate(
                                revcomp_translation_table
                            )
                    else:
                        for i in seq[2]:
                            ref_sequence += i[2]
                ref_counter += len(seq[0])
                continue
        # find location of start codon in annotated_seq v. reference
        if not ATGs and self._find_new_start:
            if not return_protein:
                return {}
            else:
                return {}, ""
        coordinates = []
        # Frame shifts: [genomic start coordinate, genomic end coordinate, CDS-
        # level start coordinate, CDS-level end coordinate, mutation info
        # associated with frame shift]
        frame_shifts = []
        counter, ref_counter = 0, 0  # hold edited transcript level coordinates
        if self._find_new_start:
            reading_frame, start_codon, ref_atg, start_warnings = self._atg_choice(
                ATGs,
                only_novel_upstream=only_novel_upstream,
                only_downstream=only_downstream,
                only_reference=only_reference,
            )
        else:
            reading_frame = 0
            start_codon = [coding_start, ref_start, True, False, False]
            ref_atg = [coding_start, ref_start, True, False, False]
            start_warnings = ["no_annotated_start_codon"]
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
        if len(start_warnings) > 0:
            transcript_warnings.append(";".join(start_warnings))
        if reading_frame:
            metadata = []
            for mutation in start_codon[2]:
                for metadata_set in mutation[2]:
                    metadata.append(metadata_set)
            frame_shifts.append([start, -1, 0, -1, metadata])
        ref_start = start_codon[1]
        coding_start = start_codon[0]
        annotated_seq.pop()
        for seq in annotated_seq:
            # skip sequence fragments that are not to be reported
            if seq[1] == "R":
                counter += len(seq[0])
                ref_counter += len(seq[0])
                continue
            elif (seq[1] == "S" and include_somatic == 2) or (
                seq[1] == "G" and include_germline == 2
            ):
                if seq[2][0][4] == "V":
                    counter += len(seq[0])
                    ref_counter += len(seq[0])
                elif seq[2][0][4] == "I":
                    counter += len(seq[0])
                    ref_counter += len(seq[0])
                #                elif seq[2][0][4] == 'D':
                #                    continue
                continue
            # skip sequence fragments that occur prior to start codon
            # handle cases where variant involves start codon
            if counter < coding_start + 2:
                if seq[1] == "H":
                    if len(seq[2][0][3]) + counter >= coding_start:
                        coordinates.append(
                            [
                                start,
                                seq[2][0][1] + len(seq[2][0][3]) * strand - 1,
                                counter,
                                counter + len(seq[2][0][3]) - 1,
                                "NA",
                                "NA",
                                seq[2][0][4],
                            ]
                        )
                    counter += len(seq[2][0][3])
                    ref_counter += len(seq[2][1][3])

                    continue
                elif seq[2][0][4] == "D":
                    ref_counter += len(seq[2][0][2])
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                start,
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
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                start,
                                seq[3] + len(seq[0]) * strand - 1,
                                0,
                                counter + len(seq[0]) - coding_start - 1,
                                "NA",
                                "NA",
                                seq[2],
                            ]
                        )
                    counter += len(seq[0])
                    if (seq[1] == "G" and include_germline == 2) or (
                        seq[1] == "S" and include_somatic == 2
                    ):
                        ref_counter += len(seq[0])
                    continue
                elif seq[2][0][4] == "V":
                    if counter + len(seq[0]) >= coding_start:
                        coordinates.append(
                            [
                                start,
                                seq[3] + len(seq[0]) * strand - 1,
                                0,
                                counter + len(seq[0]) - coding_start - 1,
                                0,
                                ref_counter + len(seq[0]) - ref_start - 1,
                                seq[2],
                            ]
                        )
                    counter += len(seq[0])
                    ref_counter += len(seq[0])
                    continue
                else:
                    # other variant types not handled at this time
                    break
            # log variants
            # handle potential frame shifts from indels
            if seq[1] == "H":
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
                    break
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
                                break
                        reading_frame = 0
                    else:
                        frame_shifts.append(
                            [seq[2][0][1], -1, counter, -1, seq[2][0][4]]
                        )
                        reading_frame = (reading_frame + read_frame1 - read_frame2) % 3
                continue
            elif seq[2][0][4] == "D":
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
                    break
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
                                break
                        reading_frame = 0
                    else:
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                        reading_frame = (reading_frame + read_frame1 - read_frame2) % 3
                ref_counter += len(seq[2][0][2])
            # handle potential frame shifts from insertions
            elif seq[2][0][4] == "I":
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
                                break
                        reading_frame = 0
                    else:
                        frame_shifts.append([seq[2][0][1], -1, counter, -1, seq[2]])
                        reading_frame = (reading_frame + len(seq[0])) % 3
                counter += len(seq[0])
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
                    A = seq_to_peptide(sequence[(i + A1) : (i + A1 + 3)])
                    B = seq_to_peptide(ref_sequence[(i + A2) : (i + A2 + 3)])
                    if A != B:
                        if frame_shifts == []:
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
                counter += len(seq[0])
                ref_counter += len(seq[0])
        # frame shifts (if they exist) continue to end of transcript
        if reading_frame != 0:
            for i in range(len(frame_shifts), 0, -1):
                if frame_shifts[i - 1][1] < 0:
                    frame_shifts[i - 1][1] = seq[3] + len(seq[0])
                    frame_shifts[i - 1][3] = counter
                else:
                    break
        protein = seq_to_peptide(sequence[start_codon[0] :], reverse_strand=False)
        protein_ref = seq_to_peptide(ref_sequence[ref_atg[1] :], reverse_strand=False)
        if '?' in protein or '?' in protein_ref:
            unknown_aa = True
        if TAA_TGA_TAG == []:
            if "X" in protein:
                for i in range(coding_start, len(sequence), 3):
                    if sequence[i : i + 3] in ["TAA", "TGA", "TAG"]:
                        coding_stop = i + 3
            else:
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
                coding_stop = len(sequence) - len(sequence) % 3
                transcript_warnings.append("nonstop")
        if len(protein) > (coding_stop - coding_start) // 3:
            if TAA_TGA_TAG != []:
                frame_shifts.append(
                    [
                        None,
                        None,
                        coding_stop,
                        3 * len(protein) + coding_start,
                        TAA_TGA_TAG[2],
                    ]
                )
        peptide_seqs = collections.defaultdict(list)
        if transcript_warnings == []:
            transcript_warnings = ("NA",)
        else:
            transcript_warnings = (";".join(transcript_warnings),)
        # get amino acid ranges for kmerization
        for size in range(min_size, max_size + 1):
            epitope_coords = []
            peptides_ref = kmerize_peptide(protein_ref, min_size=size, max_size=size)
            for coords in coordinates:
                if coords[4] != "NA":
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
                    protein[coords[0] : coords[1]], min_size=size, max_size=size
                )
                if coords[2] != "NA":
                    paired_peptides = kmerize_peptide(
                        protein_ref[coords[2] : coords[3]], min_size=size, max_size=size
                    )
                    peptide_pairs = zip(peptides, paired_peptides)
                    for pair in peptide_pairs:
                        if pair[0] not in peptides_ref:
                            if len(coords[4]) == 2 and type(coords[4][0]) == list:
                                # Dealing with peptide resulting from hybrid interval
                                data_set = coords[4][0][4]
                            else:
                                # Dealing with regular peptide
                                data_set = coords[4]
                            for mutation_data in data_set:
                                if unknown_aa and '?' in pair[0] or '?' in pair[1]:
                                    mutation_data = (
                                        mutation_data + (pair[1],) + (';'.join([transcript_warnings[0], 
                                                                                'unknown_amino_acid']),)
                                )
                                else:
                                    mutation_data = (
                                        mutation_data + (pair[1],) + transcript_warnings
                                    )
                                peptide_seqs[pair[0]].append(mutation_data)
                            peptide_seqs[pair[0]] = list(set(peptide_seqs[pair[0]]))
                else:
                    peptides = list(set(peptides).difference(peptides_ref))
                    for pep in peptides:
                        if len(coords[4]) == 2 and type(coords[4][0]) == list:
                            # Dealing with peptide resulting from hybrid interval
                            data_set = coords[4][0][4]
                        else:
                            # Dealing with regular peptide
                            data_set = coords[4]
                        for mutation_data in data_set:
                            if unknown_aa and '?' in pep:
                                mutation_data = (
                                    mutation_data + ("NA",) + (';'.join([transcript_warnings[0], 
                                                                         'unknown_amino_acid']),)
                                )
                            else:
                                mutation_data = (
                                    mutation_data + ("NA",) + transcript_warnings
                                )
                            peptide_seqs[pep].append(mutation_data)
                        peptide_seqs[pep] = list(set(peptide_seqs[pep]))
        if not return_protein:
            # return list of unique neoepitope sequences
            return peptide_seqs
        else:
            # return list of unique neoepitope sequences plus whole protein
            return peptide_seqs, protein


def gtf_to_cds(gtf_file, dictdir, pickle_it=True):
    """ References cds_dict to get cds bounds for later Bowtie query
        Keys in the dictionary are transcript IDs, while entries are lists of
            relevant CDS/stop codon data
            Data: [chromosome, sequence type, start, stop, 
                    +/- strand, transcript type]
        Writes cds_dict as a pickled dictionary
        gtf_file: input gtf file to process
        dictdir: path to directory to store pickled dicts
        Return value: dictionary
    """
    cds_dict = collections.defaultdict(list)
    cds_lines = collections.defaultdict(list)
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
                    cds_lines[transcript_id].append(tokens)
    # Sort cds_dict coordinates (left -> right) for each transcript
    for transcript_id in cds_dict:
        current_cds = cds_lines[transcript_id]
        cds_dict[transcript_id].sort(key=lambda x: x[0])
        seq_types = [x[1] for x in cds_dict[transcript_id]]
        if "start_codon" not in seq_types:
            # Fake a start codon if we have strand info
            try:
                reverse_strand = current_cds[0][6] == "-"
            except IndexError:
                # Remove incompletely annotated transcript
                del cds_dict[transcript_id]
            else:
                if reverse_strand:
                    current_cds.sort(key=lambda x: int(x[4]), reverse=True)
                    pos = int(current_cds[0][4]) - int(current_cds[0][7])
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
    # Write to pickled dictionary
    if pickle_it:
        pickle_dict = os.path.join(dictdir, "transcript_to_CDS.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(cds_dict, f)
    return cds_dict


def cds_to_tree(cds_dict, dictdir, pickle_it=True):
    """ Creates searchable tree of chromosome intervals from CDS dictionary
        Each chromosome is stored in the dictionary as an interval tree object
            Intervals are added for each CDS, with the associated transcript ID
            Assumes transcript is all on one chromosome - does not work for
                gene fusions
        Writes the searchable tree as a pickled dictionary
        cds_dict: CDS dictionary produced by gtf_to_cds()
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
    if pickle_it:
        pickle_dict = os.path.join(dictdir, "intervals_to_transcript.pickle")
        with open(pickle_dict, "wb") as f:
            pickle.dump(searchable_tree, f)
    return searchable_tree


def get_transcripts_from_tree(chrom, start, stop, cds_tree):
    """ Uses cds tree to btain transcript IDs from genomic coordinates

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
    cds = list(cds_tree[chrom].search(start, stop))
    for cd in cds:
        transcript_ids.add(cd.data)
    return list(transcript_ids)


def process_haplotypes(hapcut_output, interval_dict, phasing):
    """ Stores all haplotypes relevant to different transcripts as a dictionary
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
    try:
        if hapcut_output == "-":
            input_stream = sys.stdin
        else:
            input_stream = open(hapcut_output)
        block_transcripts = collections.defaultdict(list)
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
                        for mut in block_transcripts[transcript_id]:
                            affected_transcripts[transcript_id].append([mut])
                # Reset transcript dictionary
                block_transcripts = collections.defaultdict(list)
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
                    if alternatives[i] == "<DEL>":
                        mutation_type = "D"
                        deletion_size = len(tokens[5])
                        ref = tokens[5]
                        alt = deletion_size
                        end = pos + deletion_size
                    elif len(tokens[5]) == len(alternatives[i]):
                        mutation_type = "V"
                        pos = int(tokens[4])
                        ref = tokens[5]
                        alt = alternatives[i]
                        mut_size = len(tokens[5])
                        end = pos + mut_size
                    elif len(tokens[5]) > len(alternatives[i]):
                        mutation_type = "D"
                        deletion_size = len(tokens[5]) - len(alternatives[i])
                        pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                        ref = tokens[5][len(alternatives[i]):]
                        alt = deletion_size
                        end = pos + deletion_size
                    elif len(tokens[5]) < len(alternatives[i]):
                        mutation_type = "I"
                        insertion_size = len(alternatives[i]) - len(tokens[5])
                        pos = int(tokens[4]) + len(tokens[5]) - 1
                        ref = ""
                        alt = alternatives[i][len(tokens[5]) :]
                        end = pos + 1
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
    finally:
        if input_stream is not sys.stdin:
            input_stream.close()
    return affected_transcripts


def get_peptides_from_transcripts(
    relevant_transcripts,
    vaf_pos,
    cds_dict,
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
    include_germline=2,
    include_somatic=1,
    protein_fasta=False,
):
    """ For transcripts that are affected by a mutation, mutations are applied
        and neoepitopes resulting from mutations are called

        relevant_transcripts: dictionary linking haplotypes to transcripts;
            output from process_haplotypes()
        vaf_pos: position of VAF in VCF mutation data from HapCUT2
        cds_dict: dictionary linking transcript IDs, to lists of
            relevant CDS/stop codon data; output from gtf_to_cds()
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
        nmd: whether to include nonsense mediated decay transcripts (boolean)
        pp: whether to include polymorphic pseudogene transcripts (boolean)
        igv: whether to include IGV transcripts (boolean)
        trv: whether to include TRV transcripts (boolean)
        return value: dictionary linking neoepitopes to their associated
            metadata
        """
    neoepitopes = collections.defaultdict(list)
    fasta_entries = collections.defaultdict(set)
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
        transcript_a = Transcript(
            reference_index,
            [
                [str(chrom), "blah", seq_type, str(start), str(end), ".", strand]
                for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[
                    affected_transcript
                ]
            ],
            affected_transcript,
        )
        transcript_b = Transcript(
            reference_index,
            [
                [str(chrom), "blah", seq_type, str(start), str(end), ".", strand]
                for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[
                    affected_transcript
                ]
            ],
            affected_transcript,
        )
        # Iterate over haplotypes associated with this transcript
        haplotypes = relevant_transcripts[affected_transcript]
        for ht in haplotypes:
            # Make edits for each mutation
            for mutation in ht:
                # Determine if mutation is somatic or germline
                if mutation[6][-1] == "*":
                    mutation_class = "G"
                else:
                    mutation_class = "S"
                # Determine VAF if available
                if vaf_pos is not None:
                    try:
                        vaf_entry = mutation[6].strip("*").split(":")[vaf_pos]
                        if "," in vaf_entry:
                            vaf_entry = [x for x in vaf_entry.split(",") if x != "."]
                            if len(vaf_entry) > 0:
                                vaf = sum(
                                    [float(x.strip("%")) for x in vaf_entry]
                                ) / len(vaf_entry)
                            else:
                                vaf = None
                        else:
                            if vaf_entry.strip("%") != ".":
                                vaf = float(vaf_entry.strip("%"))
                            else:
                                vaf = None
                    except IndexError:
                        vaf = None
                else:
                    vaf = None
                # Determine which copies variant exists on & make edits
                if mutation[4] == "1":
                    transcript_a.edit(
                        mutation[3],
                        mutation[1],
                        mutation_type=mutation[7],
                        mutation_class=mutation_class,
                        vaf=vaf,
                    )
                if mutation[5] == "1":
                    transcript_b.edit(
                        mutation[3],
                        mutation[1],
                        mutation_type=mutation[7],
                        mutation_class=mutation_class,
                        vaf=vaf,
                    )
            # Extract neoepitopes
            peptides_a, protein_a = transcript_a.neopeptides(
                min_size=size_list[0],
                max_size=size_list[-1],
                include_somatic=include_somatic,
                include_germline=include_germline,
                only_novel_upstream=only_novel_upstream,
                only_downstream=only_downstream,
                only_reference=only_reference,
                return_protein=True,
            )
            peptides_b, protein_b = transcript_b.neopeptides(
                min_size=size_list[0],
                max_size=size_list[-1],
                include_somatic=include_somatic,
                include_germline=include_germline,
                only_novel_upstream=only_novel_upstream,
                only_downstream=only_downstream,
                only_reference=only_reference,
                return_protein=True,
            )
            # Store neoepitopes and their metadata
            for pep in peptides_a:
                for meta_data in peptides_a[pep]:
                    adj_meta_data = meta_data + (transcript_a.transcript_id,)
                    if adj_meta_data not in neoepitopes[pep]:
                        neoepitopes[pep].append(adj_meta_data)
            for pep in peptides_b:
                for meta_data in peptides_b[pep]:
                    adj_meta_data = meta_data + (transcript_b.transcript_id,)
                    if adj_meta_data not in neoepitopes[pep]:
                        neoepitopes[pep].append(adj_meta_data)
            if protein_fasta:
                if len(peptides_a) > 0 or len(peptides_b) > 0:
                    if protein_a == protein_b:
                        fasta_entries[affected_transcript].add(protein_a)
                    else:
                        if protein_a != "":
                            fasta_entries[affected_transcript].add(protein_a)
                        if protein_b != "":
                            fasta_entries[affected_transcript].add(protein_b)
            transcript_a.reset(reference=True)
            transcript_b.reset(reference=True)
    return neoepitopes, fasta_entries
