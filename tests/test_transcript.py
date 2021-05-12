#!/usr/bin/env python
# coding=utf-8
"""
test_transcript.py

Tests functions in transcript.py.

The MIT License (MIT)
Copyright (c) 2018 Mary A. Wood, Austin Nguyen,
                   Abhinav Nellore, and Reid Thompson

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
from inspect import getsourcefile
import os
import sys

from neoepiscope import *  # Import package in same directory as tests

import unittest

unittest.TestCase.maxDiff = None


class TestTranscript(unittest.TestCase):
    """Tests transcript object construction"""

    def setUp(self):
        """Sets up gtf file and creates dictionaries for tests"""
        self.gtf = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "tests",
            "Chr11.gtf",
        )
        self.gtf_hg38 = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "tests",
            "grch38_chr11.gtf",
        )
        self.cds, self.tx_data = gtf_to_cds(self.gtf)
        self.cds_hg38, self.tx_data_hg38 = gtf_to_cds(self.gtf_hg38)
        self.ref_prefix = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "tests",
            "Chr11.ref",
        )
        self.ref_prefix_hg38 = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "tests",
            "grch38_chr11",
        )
        self.atoi = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "neoepiscope",
                 "transcript_to_editing_hg38.pickle")
        self.rna_dict = pickle.load(open(self.atoi, "rb"))
        self.reference_index = bowtie_index.BowtieIndexReference(self.ref_prefix)
        self.reference_index_hg38 = bowtie_index.BowtieIndexReference(self.ref_prefix_hg38)
        ## All following transcripts from GRCh37 genome build ##
        # HBB-001: 628bp transcript w/ 3 exons (all coding) --> 147aa peptide
        self.transcript = Transcript(
            self.reference_index,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds[
                    "ENST00000335295.4_1"
                ]
            ],
            "ENST00000335295.4_1",
            False,
        )
        # PTDSS2-001: 2,445bp transcript w/ 12 exons (all coding) --> 487aa peptide
        self.fwd_transcript = Transcript(
            self.reference_index,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds[
                    "ENST00000308020.5_1"
                ]
            ],
            "ENST00000308020.5_1",
            False,
        )
        # OR52N1-001: 963bp transcript w/ 1 exon (all coding) --> 320aa peptide
        self.all_coding_transcript = Transcript(
            self.reference_index,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds[
                    "ENST00000317078.1_1"
                ]
            ],
            "ENST00000317078.1_1",
            False,
        )
        # CAPRIN1-001: 4108bp transcript w/ 19 exons (18/19 coding) --> 709aa peptide
        self.partial_coding_transcript = Transcript(
            self.reference_index,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds[
                    "ENST00000341394.8_1"
                ]
            ],
            "ENST00000341394.8_1",
            False,
        )
        # AP003733.1: transcript w/ 1 exon (coding)  
        self.atoi_transcript = Transcript(
            self.reference_index_hg38,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds_hg38[
                    "ENST00000318950.10"
                ]
            ],
            "ENST00000318950.10",
            False,
            self.rna_dict,
        )

        self.atoi_manual_transcript = Transcript(
            self.reference_index_hg38,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in self.cds_hg38[
                    "ENST00000318950.10"
                ]
            ],
            "ENST00000318950.10",
            False,
        )
        for item in [9750162, 9750152, 9750230, 9750194, 9750220, 9750151, 9750185, 9750221]:
            self.atoi_manual_transcript.edit("I", item, mutation_type="R", mutation_class="R", vaf=None)

        # NEAT1-002: 1745bp transcript w/ 2 exon (both non-coding) --> lncRNA
        self.non_coding_transcript = Transcript(
            self.reference_index,
            [
                [
                    str(chrom).replace("chr", ""),
                    "N/A",
                    seq_type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                ]
                for (chrom, seq_type, start, end, strand, tx_type) in [
                    ["chr11", "exon", 65190245, 65190854, "+", "lincRNA"],
                    ["chr11", "exon", 65191098, 65192298, "+", "lincRNA"],
                ]
            ],
            "ENST00000499732.2_1",
            False,
        )

    def test_transcript_structure(self):
        """Fails if structure of unedited transcript is incorrect"""
        self.assertEqual(len(self.transcript.annotated_seq()), 3)
        self.assertEqual(len(self.transcript.annotated_seq()[0][0]), 142)
        self.assertEqual(len(self.transcript.annotated_seq()[1][0]), 223)
        self.assertEqual(len(self.transcript.annotated_seq()[2][0]), 263)
        self.assertEqual(self.transcript.annotated_seq()[0][1], "R")
        self.assertEqual(
            self.transcript.intervals,
            [5246692, 5246955, 5247805, 5248028, 5248158, 5248300],
        )
        self.assertEqual(self.transcript.start_codon, 5248249)
        self.assertEqual(self.transcript._start_codon, 5248248)
        self.assertEqual(self.transcript.stop_codon, 5246828)
        self.assertEqual(self.transcript._stop_codon, 5246827)
        self.assertTrue(self.transcript.rev_strand)
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(self.transcript.deletion_intervals, [])
        self.assertEqual(self.transcript.transcript_id, "ENST00000335295.4_1")

    # Reading frame tests
    def test_rev_reading_frame(self):
        """Fails if incorrect reading frame is called in rev transcript"""
        # Before and after exons
        self.assertIsNone(self.transcript.reading_frame(5248350))
        self.assertIsNone(self.transcript.reading_frame(5246600))
        # In exons, before or after start/stop codons
        self.assertIsNone(self.transcript.reading_frame(5248290))
        self.assertIsNone(self.transcript.reading_frame(5246817))
        # In start codon exon
        self.assertEqual(self.transcript.reading_frame(5248161), 0)
        self.assertEqual(self.transcript.reading_frame(5248217), 1)
        self.assertEqual(self.transcript.reading_frame(5248249), 2)
        # In other exon
        self.assertEqual(self.transcript.reading_frame(5246848), 0)
        self.assertEqual(self.transcript.reading_frame(5247859), 1)
        self.assertEqual(self.transcript.reading_frame(5248014), 2)

    def test_fwd_reading_frame(self):
        """Fails if incorrect reading frame is called in fwd transcript"""
        # Before and after exons
        self.assertIsNone(self.fwd_transcript.reading_frame(450200))
        self.assertIsNone(self.fwd_transcript.reading_frame(491392))
        # In exons, before or after start/stop codons
        self.assertIsNone(self.fwd_transcript.reading_frame(450394))
        self.assertIsNone(self.fwd_transcript.reading_frame(490824))
        # In start codon exon
        self.assertEqual(self.fwd_transcript.reading_frame(450459), 0)
        self.assertEqual(self.fwd_transcript.reading_frame(450550), 1)
        self.assertEqual(self.fwd_transcript.reading_frame(450578), 2)
        # In other exon
        self.assertEqual(self.fwd_transcript.reading_frame(487029), 0)
        self.assertEqual(self.fwd_transcript.reading_frame(460270), 1)
        self.assertEqual(self.fwd_transcript.reading_frame(460259), 2)
        # In non-coding transcript
        self.assertIsNone(self.non_coding_transcript.reading_frame(65190449))
        # In a transcript with one or more non-coding exons
        self.assertIsNone(self.partial_coding_transcript.reading_frame(34073265))
        self.assertEqual(self.partial_coding_transcript.reading_frame(34073973), 2)

    # Transcript sequence editing tests
    def test_irrelevant_edit(self):
        """Fails if edit is made for intronic position"""
        # In transcript with all coding exons
        self.transcript.edit("A", 5248155)
        relevant_edits = self.transcript.expressed_edits()
        self.assertEqual(
            self.transcript.edits[5248154],
            [("A", "V", "S", ("11", 5248155, "C", "A", "V", None))],
        )
        self.assertEqual(relevant_edits[0], {})
        self.assertEqual(
            relevant_edits[1],
            [
                (5246692, "R", ()),
                (5246955, "R", ()),
                (5247805, "R", ()),
                (5248028, "R", ()),
                (5248158, "R", ()),
                (5248300, "R", ()),
            ],
        )
        self.assertEqual(len(self.transcript.annotated_seq()), 3)
        self.assertEqual(
            len([x for x in self.transcript.annotated_seq() if x[1] != "R"]), 0
        )
        # Insertion in an intron, only fetching a portion of annotated seq
        self.transcript.reset()
        self.assertEqual(
            self.transcript.annotated_seq(5246947, 5247012),
            [("CTCCTGGGCA", "R", [()], 5246956)],
        )
        self.transcript.edit("ATGCGGGG", 5246980, mutation_type="I")
        self.assertEqual(
            self.transcript.annotated_seq(5246947, 5247012),
            [("CTCCTGGGCA", "R", [()], 5246956)],
        )
        # In a transcript with no coding exons
        self.non_coding_transcript.edit("A", 65190856)
        relevant_edits = self.non_coding_transcript.expressed_edits()
        self.assertEqual(
            self.non_coding_transcript.edits[65190855],
            [("A", "V", "S", ("11", 65190856, "T", "A", "V", None))],
        )
        self.assertEqual(relevant_edits[0], {})
        self.assertEqual(
            relevant_edits[1],
            [
                (65190243, "R", ()),
                (65190853, "R", ()),
                (65191096, "R", ()),
                (65192297, "R", ()),
            ],
        )
        # In a coding transcript with a noncoding exon
        self.partial_coding_transcript.edit("T", 34073965)
        relevant_edits = self.partial_coding_transcript.expressed_edits()
        self.assertEqual(
            self.partial_coding_transcript.edits[34073964],
            [("T", "V", "S", ("11", 34073965, "A", "T", "V", None))],
        )
        self.assertEqual(relevant_edits[0], {})
        self.assertEqual(len(relevant_edits[1]), 38)

    def test_relevant_edit(self):
        """Fails if edit is not made for position within exon"""
        # In transcript with all coding exons
        self.transcript.edit("A", 5248299)
        relevant_edits = self.transcript.expressed_edits()
        self.assertEqual(
            self.transcript.edits[5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(
            relevant_edits[0][5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(
            relevant_edits[1],
            [
                (5246692, "R", ()),
                (5246955, "R", ()),
                (5247805, "R", ()),
                (5248028, "R", ()),
                (5248158, "R", ()),
                (5248300, "R", ()),
            ],
        )
        self.assertEqual(len(self.transcript.annotated_seq()), 5)
        self.assertEqual(self.transcript.annotated_seq()[1][0], "T")
        self.assertEqual(self.transcript.annotated_seq()[1][1], "S")
        self.assertEqual(
            len([x for x in self.transcript.annotated_seq() if x[1] != "R"]), 1
        )
        # In a transcript with no coding exons
        self.non_coding_transcript.edit("A", 65190277)
        relevant_edits = self.non_coding_transcript.expressed_edits()
        self.assertEqual(
            relevant_edits[0],
            {65190276: [("A", "V", "S", ("11", 65190277, "C", "A", "V", None))]},
        )
        # In a coding transcript with a noncoding exon
        self.partial_coding_transcript.edit("A", 34073237)
        relevant_edits = self.partial_coding_transcript.expressed_edits()
        self.assertEqual(
            relevant_edits[0],
            {34073236: [("A", "V", "S", ("11", 34073237, "C", "A", "V", None))]},
        )

    def test_reset_to_reference(self):
        """Fails if transcript is not reset to reference"""
        self.transcript.edit("A", 5248299)
        self.transcript.reset(reference=True)
        self.assertEqual(self.transcript.edits, {})

    def test_edit_and_save(self):
        """Fails if edits aren't saved"""
        self.transcript.edit("A", 5248299)
        self.transcript.edit(3, 5246694, mutation_type="D")
        self.transcript.save()
        self.assertEqual(
            self.transcript.last_edits[5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(
            self.transcript.last_deletion_intervals,
            [(5246692, 5246695, "S", ("11", 5246694, "TTG", "", "D", None, False))],
        )

    def test_reset_to_save_point(self):
        """Fails if new edit not erased or old edits not retained"""
        self.transcript.edit("A", 5248299)
        self.transcript.edit(3, 5246694, mutation_type="D")
        self.transcript.save()
        self.transcript.edit("G", 5248165)
        self.assertEqual(
            self.transcript.edits[5248164],
            [("G", "V", "S", ("11", 5248165, "C", "G", "V", None))],
        )
        self.transcript.reset(reference=False)
        self.assertNotIn(2182387, self.transcript.edits)
        self.assertEqual(
            self.transcript.last_edits[5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(
            self.transcript.last_deletion_intervals,
            [(5246692, 5246695, "S", ("11", 5246694, "TTG", "", "D", None, False))],
        )
        self.assertNotEqual(self.transcript.edits, {})

    def test_SNV_seq(self):
        """Fails if SNV is edited incorrectly"""
        self.transcript.edit("A", 5248299)
        self.assertEqual(
            self.transcript.edits[5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(seq[0], ("AC", "R", [()], 5248301))
        self.assertEqual(
            seq[1], ("T", "S", [("11", 5248299, "T", "A", "V", None)], 5248299)
        )
        self.assertEqual(len(seq[2][0]), 139)
        self.assertEqual(len(seq[3][0]), 223)
        self.assertEqual(len(seq[4][0]), 263)

    def test_inside_insertion(self):
        """Fails if indel within one exon is inserted incorrectly"""
        self.transcript.edit("Q", 5248165, mutation_type="I")
        self.assertEqual(
            self.transcript.edits[5248164],
            [("Q", "I", "S", ("11", 5248165, "", "Q", "I", None))],
        )
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 6)
        self.assertEqual(
            seq[1], ("Q", "S", [("11", 5248165, "", "Q", "I", None)], 5248165)
        )
        self.assertEqual(len(seq[0][0]), 136)
        self.assertEqual(len(seq[2][0]), 1)
        self.assertEqual(len(seq[3][0]), 5)
        self.assertEqual(len(seq[4][0]), 223)
        self.assertEqual(len(seq[5][0]), 263)

    def test_adjacent_indel(self):
        """Fails if exon-adjacent indel is inserted incorrectly"""
        self.transcript.edit("Q", 5248029, mutation_type="I")
        self.assertEqual(
            self.transcript.edits[5248028],
            [("Q", "I", "S", ("11", 5248029, "", "Q", "I", None))],
        )
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(
            seq[1], ("Q", "S", [("11", 5248029, "", "Q", "I", None)], 5248029)
        )
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[2][0]), 1)
        self.assertEqual(len(seq[3][0]), 222)
        self.assertEqual(len(seq[4][0]), 263)

    def test_inside_deletion(self):
        """Fails if deletion completely within an exon is made improperly"""
        self.transcript.edit(3, 5246700, mutation_type="D")
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(
            self.transcript.deletion_intervals,
            [(5246698, 5246701, "S", ("11", 5246700, "TGA", "", "D", None, False))],
        )
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(
            seq[3], ("", "S", [("11", 5246700, "TGA", "", "D", None)], 5246700)
        )
        self.assertEqual(seq[4], ("TTGCAA", "R", [()], 5246699))
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[1][0]), 223)
        self.assertEqual(len(seq[2][0]), 254)

    def test_overlapping_deletion(self):
        """Fails if deletion of an exon-intron junction is incorrect"""
        self.transcript.edit(10, 5246950, mutation_type="D")
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(
            self.transcript.deletion_intervals,
            [
                (
                    5246948,
                    5246958,
                    "S",
                    ("11", 5246950, "CCAGGAGCTG", "", "D", None, True),
                )
            ],
        )
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 2)
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[1][0]), 223)

    def test_complex_overlapping_deletion(self):
        """Fails if deletion of an exon-intron junction + other mutations is incorrect"""
        self.transcript.edit(10, 5247804, mutation_type="D")
        self.transcript.edit("G", 5247943)
        self.assertEqual(
            self.transcript.deletion_intervals,
            [
                (
                    5247802,
                    5247812,
                    "S",
                    ("11", 5247804, "CACCCTGAAG", "", "D", None, True),
                )
            ],
        )
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 1)
        self.assertEqual(len(seq[0][0]), 142)

    def test_spanning_deletion(self):
        """Fails if deletion spanning two or more exons is incorrect"""
        # Spanning 2 exons
        self.transcript.edit(137, 5248025, mutation_type="D")
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(
            self.transcript.deletion_intervals[0][0:3], (5248023, 5248160, "S")
        )
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 4)
        self.assertEqual(len(seq[1][2][0][2]), 137)
        self.assertEqual(seq[1][0:2], ("", "S"))
        self.assertEqual(len(seq[0][0]), 140)
        self.assertEqual(len(seq[2][0]), 218)
        self.assertEqual(len(seq[3][0]), 263)
        # Spanning more than 2 exons
        self.partial_coding_transcript.edit(23935, 34074178, mutation_type="D")
        seq = self.partial_coding_transcript.annotated_seq()
        self.assertEqual(len(seq), 17)
        self.assertEqual(len(seq[2][2][0][2]), 23935)
        self.assertEqual(len(seq[1][0]), 210)
        self.assertEqual(len(seq[3][0]), 77)

    def test_deletion_over_transcript_start(self):
        """Fails if deletion spanning start of transcript is incorrect"""
        self.transcript.edit(5, 5246692, mutation_type="D")
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 2)
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[1][0]), 223)
        self.fwd_transcript.edit(10, 450275, mutation_type="D")
        seq2 = self.fwd_transcript.annotated_seq()
        self.assertEqual(seq2, [])

    def test_adjacent_deletions(self):
        """Fails if adjacent deletions are handled incorrectly"""
        self.fwd_transcript.edit(5, 450286, mutation_type="D", mutation_class="G")
        self.fwd_transcript.edit(5, 450291, mutation_type="D")
        seq = self.fwd_transcript.annotated_seq()
        self.assertEqual(len(seq), 15)
        self.assertEqual(seq[1][1], "G")
        self.assertEqual(seq[2][1], "S")
        self.transcript.edit(5, 5248292, mutation_type="D", mutation_class="G")
        self.transcript.edit(5, 5248287, mutation_type="D")
        seq2 = self.transcript.annotated_seq()
        self.assertEqual(len(seq2), 6)
        self.assertEqual(seq[1][1], "G")
        self.assertEqual(seq[2][1], "S")

    def test_hybrid_deletions(self):
        """Fails if overlapping germline and somatic deletions are hybridized
        incorrectly"""
        # Forward transcript
        self.fwd_transcript.edit(5, 450550, mutation_type="D")
        self.fwd_transcript.edit(5, 450552, mutation_type="D", mutation_class="G")
        seq = self.fwd_transcript.annotated_seq()
        self.assertEqual(seq[1][0], "")
        self.assertEqual(seq[1][1], "H")
        self.assertEqual(
            seq[1][2][0],
            [
                "11",
                450550,
                "CGTCTGC",
                "",
                [
                    ("11", 450550, "CGTCT", "", "D", None),
                    ("11", 450552, "TCTGC", "", "D", None),
                ],
            ],
        )
        self.assertEqual(
            seq[1][2][1],
            ["11", 450552, "TCTGC", "CG", [("11", 450552, "TCTGC", "", "D", None)]],
        )
        # Reverse transcript
        self.transcript.edit(5, 5248211, mutation_type="D", mutation_class="G")
        self.transcript.edit(5, 5248208, mutation_type="D")
        seq2 = self.transcript.annotated_seq()
        self.assertEqual(seq2[1][0], "")
        self.assertEqual(seq2[1][1], "H")
        self.assertEqual(
            seq2[1][2][0],
            [
                "11",
                5248208,
                "AGGGCAGT",
                "",
                [
                    ("11", 5248211, "GCAGT", "", "D", None),
                    ("11", 5248208, "AGGGC", "", "D", None),
                ],
            ],
        )
        self.assertEqual(
            seq2[1][2][1],
            ["11", 5248211, "GCAGT", "AGG", [("11", 5248211, "GCAGT", "", "D", None)]],
        )

    def test_hybrid_boundary_deletion(self):
        """Fails if hybrid deletion over intron-exon boundary is incorrect"""
        self.fwd_transcript.edit(4, 473888, mutation_type="D", mutation_class="G")
        self.fwd_transcript.edit(10, 473890, mutation_type="D")
        seq = self.fwd_transcript.annotated_seq()
        self.assertEqual(len(seq), 2)
        self.transcript.edit(10, 5246952, mutation_type="D", mutation_class="G")
        self.transcript.edit(20, 5246955, mutation_type="D")
        seq2 = self.transcript.annotated_seq()
        self.assertEqual(len(seq2), 2)

    def test_complex_indel(self):
        """Fails if adjacent deletion and insertion are handled incorrectly"""
        self.transcript.edit(4, 5247824, mutation_type="D")
        self.transcript.edit("TT", 5247827, mutation_type="I")
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 7)
        self.fwd_transcript.edit(4, 473916, mutation_type="D")
        self.fwd_transcript.edit("AA", 473919, mutation_type="I")
        seq = self.fwd_transcript.annotated_seq()
        self.assertEqual(len(seq), 16)

    def test_deletion_of_transcript(self):
        """Fails if deletion of entire transcript is incorrect"""
        self.transcript.edit(1700, 5246690, mutation_type="D")
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 1)
        self.assertEqual(seq[0][0], "")
        self.assertEqual(len(seq[0][2][0][2]), 1700)

    def test_compound_variants(self):
        """Fails if transcript with multiple variant types is incorrect"""
        self.transcript.edit(137, 5248025, mutation_type="D")
        self.transcript.edit("Q", 5248165, mutation_type="I")
        self.transcript.edit("A", 5248299)
        self.assertEqual(len(self.transcript.edits.keys()), 2)
        self.assertEqual(
            self.transcript.edits[5248298],
            [("A", "V", "S", ("11", 5248299, "T", "A", "V", None))],
        )
        self.assertEqual(
            self.transcript.edits[5248164],
            [("Q", "I", "S", ("11", 5248165, "", "Q", "I", None))],
        )
        self.assertEqual(
            self.transcript.deletion_intervals[0][0:3], (5248023, 5248160, "S")
        )
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 9)
        self.assertEqual(seq[0], ("AC", "R", [()], 5248301))
        self.assertEqual(
            seq[1], ("T", "S", [("11", 5248299, "T", "A", "V", None)], 5248299)
        )
        self.assertEqual(len(seq[2][0]), 133)
        self.assertEqual(
            seq[3], ("Q", "S", [("11", 5248165, "", "Q", "I", None)], 5248165)
        )
        self.assertEqual(seq[4], ("G", "R", [()], 5248165))
        self.assertEqual(seq[5], ("GGC", "R", [()], 5248164))
        self.assertEqual(len(seq[6][2][0][2]), 137)
        self.assertEqual(len(seq[7][0]), 218)
        self.assertEqual(len(seq[8][0]), 263)

    # Neopeptide tests
    def test_no_mutations_peptides(self):
        """Fails if peptides are returned for unmutated sequence"""
        peptides = self.fwd_transcript.neopeptides()
        self.assertFalse(peptides)
        rev_peptides = self.transcript.neopeptides()
        self.assertFalse(rev_peptides)

    def test_noncoding_mutation_peptides(self):
        """Fails if peptides are returned for mutation in noncoding
        sequence"""
        self.fwd_transcript.edit("G", 450286)
        peptides = self.fwd_transcript.neopeptides()
        self.assertFalse(peptides)
        self.transcript.edit("A", 5248266)
        rev_peptides = self.transcript.neopeptides()
        self.assertFalse(rev_peptides)

    def test_synonymous_snv_peptides(self):
        """Fails if peptides are returned for a synonymous snv"""
        self.fwd_transcript.edit("A", 450464)
        peptides = self.fwd_transcript.neopeptides()
        self.assertFalse(peptides)
        self.transcript.edit("A", 5248005)
        rev_peptides = self.transcript.neopeptides()
        self.assertFalse(rev_peptides)

    def test_missense_snv_peptides(self):
        """Fails if incorrect peptides are returned for missense SNV"""
        self.fwd_transcript.edit("T", 450502)
        peptides = self.fwd_transcript.neopeptides().keys()
        f_peptides = [pep for pep in peptides if "F" in pep]
        self.assertEqual(len(peptides), 38)
        self.assertEqual(len(peptides), len(f_peptides))
        self.assertEqual(sorted(peptides)[0], "AGGPRPEF")
        self.assertEqual(sorted(peptides)[-1], "RRDAGGPRPEF")
        self.transcript.edit("T", 5248006)
        rev_peptides = self.transcript.neopeptides().keys()
        n_peptides = [pep for pep in rev_peptides if "N" in pep]
        self.assertEqual(len(rev_peptides), 38)
        self.assertEqual(len(rev_peptides), len(n_peptides))
        self.assertEqual(sorted(rev_peptides)[0], "GRLLVVYPWN")
        self.assertEqual(sorted(rev_peptides)[-1], "YPWNQRFFESF")

    def test_in_frame_insertion_peptides(self):
        """Fails if incorrect peptides are returned for in-frame
        insertion"""
        self.fwd_transcript.edit("AAA", 450551, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides().keys()
        k_peptides = [pep for pep in peptides if "K" in pep]
        self.assertEqual(len(peptides), 38)
        self.assertEqual(len(peptides), len(k_peptides))
        self.assertEqual(sorted(peptides)[0], "ASLEEPPDGPK")
        self.assertEqual(sorted(peptides)[-1], "SLEEPPDGPKS")
        self.transcript.edit("TTT", 5247986, mutation_type="I")
        rev_peptides = self.transcript.neopeptides().keys()
        k_peptides = [pep for pep in peptides if "K" in pep]
        self.assertEqual(len(rev_peptides), 38)
        self.assertEqual(len(rev_peptides), len(k_peptides))
        self.assertEqual(sorted(rev_peptides)[0], "ESKFGDLS")
        self.assertEqual(sorted(rev_peptides)[-1], "YPWTQRFFESK")

    def test_synonymous_inframe_insertion_peptides(self):
        """Fails if incorrect peptides are returned for an insertion into
        a codon that maintains the AA sequence of that codon"""
        self.fwd_transcript.edit("AAA", 450502, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 38)
        self.transcript.edit("TTT", 5246874, mutation_type="I")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 34)

    def test_frameshift_insertion(self):
        """Fails if incorrect peptides are returned for frameshift
        insertion"""
        self.fwd_transcript.edit("AAAAA", 473925, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 108)
        self.assertEqual(sorted(peptides)[0], "ASILVFLK")
        self.assertEqual(sorted(peptides)[-1], "VLESHKLKTGH")
        self.transcript.edit("AAAAA", 5247883, mutation_type="I")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(sorted(rev_peptides)[0], "AFSDGLAHLV")
        self.assertEqual(sorted(rev_peptides)[-1], "VFTTSRAPLPH")

    def test_in_frame_deletion_peptides(self):
        """Fails if incorrect peptides are given for in-frame deletion"""
        self.fwd_transcript.edit(3, 450555, mutation_type="D")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 34)
        self.assertEqual(sorted(peptides)[0], "DGPSGQAT")
        self.assertEqual(sorted(peptides)[-1], "SLEEPPDGPSG")
        self.transcript.edit(3, 5247858, mutation_type="D")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 34)
        self.assertEqual(sorted(rev_peptides)[0], "ALSELHCD")
        self.assertEqual(sorted(rev_peptides)[-1], "TFALSELHCDK")

    def test_synonymous_inframe_deletion_peptides(self):
        """Fails if incorrect peptides are returned for a deletion of
        a codon that maintains the AA sequence of that codon"""
        self.fwd_transcript.edit(3, 473918, mutation_type="D")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 34)
        self.transcript.edit(3, 5247922, mutation_type="D")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 30)

    def test_frameshift_deletion(self):
        """Fails if incorrect peptides are returned for frameshift
        deletion"""
        self.fwd_transcript.edit(5, 473912, mutation_type="D")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 40)
        self.assertEqual(sorted(peptides)[0], "ASSFLMFW")
        self.assertEqual(sorted(peptides)[-1], "YNTKRGIVASS")
        self.transcript.edit(5, 5247930, mutation_type="D")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 32)
        self.assertEqual(sorted(rev_peptides)[0], "AVMGNPKVKG")
        self.assertEqual(sorted(rev_peptides)[-1], "VMGNPKVKGQE")

    def test_nonstop_mutation_peptides(self):
        """Fails if mutation altering stop codon does not return peptides
        past the end of the original peptide to the new stop"""
        self.fwd_transcript.edit("A", 490580)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 768)
        self.assertEqual(sorted(peptides)[0], "AARVGARP")
        self.assertEqual(sorted(peptides)[-1], "YTGSTSTSPAA")
        self.transcript.edit("G", 5246828)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 84)
        self.assertEqual(sorted(rev_peptides)[0], "AHKYHYAR")
        self.assertEqual(sorted(rev_peptides)[-1], "YHYARFLAVQF")
        # Check warnings for transcript with stop codon removed
        self.all_coding_transcript.edit("T", 5809086)
        peptides = self.all_coding_transcript.neopeptides()
        for pep in peptides:
            for mutation_data in peptides[pep]:
                self.assertEqual(
                    mutation_data[7], "annotated_stop_codon_disrupted;nonstop"
                )

    def test_split_start(self):
        """Fails if split start codon is handled improperly"""
        self.fwd_transcript.edit("AT", 450419, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(peptides), 278)
        self.assertEqual(sorted(peptides)[0], "AAAAPSPR")
        self.assertEqual(sorted(peptides)[-1], "WRSRLTGRLPA")
        self.transcript.edit("T", 5248290, mutation_type="I")
        rev_peptides = self.transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(rev_peptides), 42)
        self.assertEqual(sorted(rev_peptides)[0], "ATSNRHHG")
        self.assertEqual(sorted(rev_peptides)[-1], "TSNRHHGASDS")
        self.partial_coding_transcript.edit("A", 34073379)
        partial_peptides = self.partial_coding_transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(partial_peptides), 122)
        self.assertEqual(sorted(partial_peptides)[0], "AHDALGHQ")
        self.assertEqual(sorted(partial_peptides)[-1], "VVRTATAVGFL")

    def test_start_lost_peptides(self):
        """Fails if mutation altering start codon does not return peptides
        from a new start codon"""
        self.fwd_transcript.edit("T", 450456)
        peptides = self.fwd_transcript.neopeptides()
        self.assertFalse(peptides)
        # Next start immediately followed by stop for fwd strand transcript
        self.transcript.edit("G", 5248251)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 122)
        self.assertEqual(sorted(rev_peptides)[0], "AGCWWSTL")
        self.assertEqual(sorted(rev_peptides)[-1], "WWSTLGPRGSL")
        self.partial_coding_transcript.edit("T", 34073968)
        partial_peptides = self.partial_coding_transcript.neopeptides()
        self.assertFalse(partial_peptides)
        # Next start is in the same reading frame
        self.all_coding_transcript.edit("T", 5810046)
        all_peptides = self.all_coding_transcript.neopeptides()
        self.assertFalse(all_peptides)
        # From next start, peptide is only 4 aa

    def test_start_lost_and_new_inframe_start(self):
        """Fails if peptides aren't returned from a new in frame start
        codon when the original is disrupted"""
        self.fwd_transcript.edit("ATG", 450446, mutation_type="I")
        self.fwd_transcript.edit("T", 450456)
        peptides = self.fwd_transcript.neopeptides(
            only_novel_upstream=True,
            only_downstream=False,
        ).keys()
        self.assertEqual(len(peptides), 20)
        self.transcript.edit("CAT", 5248263, mutation_type="I")
        self.transcript.edit("G", 5248251)
        rev_peptides = self.transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(rev_peptides), 24)
        self.partial_coding_transcript.edit("T", 34073968)
        self.partial_coding_transcript.edit("ATG", 34073397, mutation_type="I")
        partial_peptides = self.partial_coding_transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(partial_peptides), 36)

    def test_start_lost_and_new_out_of_frame_start(self):
        """Fails if peptides aren't returned from a new out of frame start
        codon when the original is disrupted"""
        self.fwd_transcript.edit("ATG", 450445, mutation_type="I")
        self.fwd_transcript.edit("T", 450456)
        peptides = self.fwd_transcript.neopeptides(
            only_novel_upstream=True,
            only_downstream=False,
        ).keys()
        self.assertEqual(len(peptides), 98)
        self.transcript.edit("CAT", 5248280, mutation_type="I")
        self.transcript.edit("G", 5248251)
        rev_peptides = self.transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(rev_peptides), 22)
        self.partial_coding_transcript.edit("T", 34073968)
        self.partial_coding_transcript.edit("ATG", 34073396, mutation_type="I")
        partial_peptides = self.partial_coding_transcript.neopeptides(
            only_novel_upstream=True, only_downstream=False
        ).keys()
        self.assertEqual(len(partial_peptides), 102)

    def test_skipping_new_start(self):
        """Fails if peptides are returned from a new upstream start codon
        when only searching downstream"""
        self.fwd_transcript.edit("ATG", 450445, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides(only_downstream=True)
        self.assertFalse(peptides)
        self.transcript.edit("CAT", 5248280, mutation_type="I")
        rev_peptides = self.transcript.neopeptides(only_downstream=True)
        self.assertFalse(rev_peptides)
        self.partial_coding_transcript.edit("ATG", 34073397, mutation_type="I")
        partial_peptides = self.partial_coding_transcript.neopeptides(
            only_downstream=True
        ).keys()
        self.assertFalse(partial_peptides)
        self.non_coding_transcript.edit("ATG", 65190565, mutation_type="I")
        noncoding_peptides = self.non_coding_transcript.neopeptides()
        self.assertEqual(noncoding_peptides, {})

    def test_translation_lost_and_restored(self):
        """Fails if germline start lost with somatic restore returns
        incorrect peptides"""
        self.fwd_transcript.edit("C", 450458, mutation_class="G")
        self.fwd_transcript.edit("G", 450458)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 1914)
        self.all_coding_transcript.edit("C", 5810042, mutation_class="G")
        self.all_coding_transcript.edit("G", 5810042)
        coding_peptides = self.all_coding_transcript.neopeptides().keys()
        self.assertEqual(len(coding_peptides), 1246)

    def test_compound_indel_peptides(self):
        """Fails if incorrect peptides are returned with complementary
        indels are introduced (i.e. frameshift then return to frame)"""
        self.fwd_transcript.edit(4, 473924, mutation_type="D")
        self.fwd_transcript.edit("AAAA", 473952, mutation_type="I")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 74)
        self.assertEqual(sorted(peptides)[0], "ASILVFFL")
        self.assertEqual(sorted(peptides)[-1], "VFFLESHKLKT")
        self.transcript.edit(4, 5247921, mutation_type="D")
        self.transcript.edit("AAAA", 5247933, mutation_type="I")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 50)

    def test_all_coding_tx(self):
        """Fails if transcript that is all coding sequence is
        handled improperly"""
        self.all_coding_transcript.edit("C", 5810036)
        peptides = self.all_coding_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 16)

    def test_hybrid_deletion_peptides(self):
        """Fails if peptides from hybrid deletion are returned incorrectly"""
        # Forward transcript
        self.fwd_transcript.edit(5, 450550, mutation_type="D")
        self.fwd_transcript.edit(5, 450552, mutation_type="D", mutation_class="G")
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 124)
        self.assertEqual(sorted(peptides)[0], "AAAAPSPR")
        self.assertEqual(sorted(peptides)[-1], "TTTAPTPSSGE")
        # Reverse transcript
        self.transcript.edit(5, 5248211, mutation_type="D", mutation_class="G")
        self.transcript.edit(5, 5248208, mutation_type="D")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 28)

    def test_germline_vs_somatic(self):
        """ "Fails if incorrect peptides are returned for different
        germline/somatic mutation inclusion decisions"""
        self.fwd_transcript.edit("T", 450502)
        self.fwd_transcript.edit("T", 450503, mutation_class="G")
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=0, include_germline=0
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=0, include_germline=1
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=0, include_germline=2
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=1, include_germline=0
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=1, include_germline=1
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=1, include_germline=2
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=2, include_germline=0
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=2, include_germline=1
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.fwd_transcript.neopeptides(
                    include_somatic=2, include_germline=2
                ).keys()
            ),
            0,
        )
        self.transcript.edit("T", 5248006)
        self.transcript.edit("C", 5248007, mutation_class="G")
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=0, include_germline=0
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=0, include_germline=1
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=0, include_germline=2
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=1, include_germline=0
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=1, include_germline=1
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=1, include_germline=2
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=2, include_germline=0
                ).keys()
            ),
            0,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=2, include_germline=1
                ).keys()
            ),
            38,
        )
        self.assertEqual(
            len(
                self.transcript.neopeptides(
                    include_somatic=2, include_germline=2
                ).keys()
            ),
            0,
        )

    def test_compound_all(self):
        """Fails if incorrect peptides are returned when multiple
        germline/somatic mutations are introduced"""
        # Forward transcript
        self.fwd_transcript.edit("AAA", 490579, mutation_type="I", mutation_class="G")
        self.fwd_transcript.edit("C", 490580, mutation_type="V")
        self.fwd_transcript.edit("A", 490582, mutation_type="I")
        self.fwd_transcript.edit("A", 490586, mutation_type="V")
        self.fwd_transcript.edit("G", 490587, mutation_type="D")
        self.fwd_transcript.edit("C", 490588, mutation_type="V")
        self.fwd_transcript.edit("A", 490593, mutation_type="V", mutation_class="G")
        self.fwd_transcript.edit("GAGGAGGAGGAG", 450536, mutation_type="I")
        self.fwd_transcript.save()
        self.fwd_transcript.edit("T", 450457, mutation_type="D")
        self.fwd_transcript.reset()
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 58)
        self.assertEqual(sorted(peptides)[5], "APTPNKRTYP")
        self.assertEqual(sorted(peptides)[-1], "VPAGRASLEEE")
        self.assertEqual(sorted(peptides)[-5], "SLEEEEEEP")
        peptides = self.fwd_transcript.neopeptides(include_germline=1).keys()
        self.assertEqual(len(peptides), 62)
        self.assertEqual(sorted(peptides)[0], "AEGEGAPTPNK")
        self.assertEqual(sorted(peptides)[32], "EGEGAPTPNKR")
        self.assertEqual(sorted(peptides)[-5], "SLEEEEEEP")
        # Reverse transcript
        self.transcript.edit(4, 5247921, mutation_type="D")
        self.transcript.edit("AAAA", 5247933, mutation_type="I")
        self.transcript.edit(3, 5248170, mutation_type="D", mutation_class="G")
        self.transcript.edit("T", 5246872, mutation_type="V")
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 88)
        self.assertEqual(sorted(rev_peptides)[0], "AAYQKMVA")
        self.assertEqual(sorted(rev_peptides)[-1], "YQKMVAGVANA")
        rev_peptides = self.transcript.neopeptides(include_germline=1).keys()
        self.assertEqual(len(rev_peptides), 122)
        self.assertEqual(sorted(rev_peptides)[13], "DEVGGALG")
        self.assertEqual(sorted(rev_peptides)[-13], "VNVDEVGGALG")
        self.assertEqual(sorted(rev_peptides)[0], "AAYQKMVA")
        # Partially coding transcript
        self.partial_coding_transcript.edit(
            "AAA", 34073967, mutation_type="I", mutation_class="G"
        )
        self.partial_coding_transcript.edit(1, 34073973, mutation_type="D")
        self.partial_coding_transcript.edit("T", 34073988, mutation_type="I")
        self.partial_coding_transcript.edit(
            9, 34074019, mutation_type="D", mutation_class="G"
        )
        partial_peptides = self.partial_coding_transcript.neopeptides().keys()
        self.assertEqual(len(partial_peptides), 28)
        # All coding transcript
        self.all_coding_transcript.edit(
            "AAA", 5810046, mutation_type="I", mutation_class="G"
        )
        self.all_coding_transcript.edit(1, 5810032, mutation_type="D")
        self.all_coding_transcript.edit("T", 5810012, mutation_type="I")
        all_peptides = self.all_coding_transcript.neopeptides()
        self.assertEqual(
            all_peptides,
            {"MSFLKAPA": [("11", 5810032, "A", "", "D", None, "NA", "NA")]},
        )

    def test_expressed_edits_with_rna_edits_from_dict(self):
        """check expressed_edit can read and generate edits using
            rna_editing_sites from dictionary"""
        self.atoi_transcript.expressed_edits(include_rna_edits=True)
        self.assertEqual(self.atoi_transcript.edits[9750161],
                [('I', 'R', 'R', ('11', 9750162, 'A', 'I', 'R', None))])
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=False)
        self.assertNotIn(9750161, edits)

    def test_edit_with_rna_edits_at_start_codon(self):
        """check whether rna editing handled correctly for start codon"""
        self.atoi_transcript.edit('I', 9664180, mutation_type="R", mutation_class="R", vaf=None)
        self.assertEqual(self.atoi_transcript.edits[9664179],
                [('I', 'R', 'R', ('11', 9664180, 'A', 'I', 'R', None))])
        self.assertEqual(self.atoi_transcript.all_transcript_warnings, ["rna_editing_may_disrupt_start_codon"])

    def test_expressed_rna_edit_with_not_A_in_ref_genome(self):
        pos = 9750164
        self.atoi_transcript.edit('I', pos, mutation_type="R", mutation_class="R", vaf=None)
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True)
        self.assertNotIn(pos-1, edits)

        self.atoi_manual_transcript.edit('I', pos, mutation_type="R", mutation_class="R", vaf=None)
        edits, _ = self.atoi_manual_transcript.expressed_edits(include_rna_edits=True)
        self.assertNotIn(pos-1, edits)

    def test_expressed_edit_with_overlapping_rna_edit_and_germline_with_germline_options(self):
        pos = 9750162
        self.atoi_transcript.edit('C', pos, mutation_type="V", mutation_class="G", vaf=None)
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True, include_germline=0)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True, include_germline=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "C")
        self.assertEqual(edit[3][3], "C")
        self.assertEqual(edit[3][4], "V")


    def test_expressed_edit_with_overlapping_rna_edit_and_somatic_with_somatic_options(self):
        pos = 9750162
        self.atoi_transcript.edit('C', pos, mutation_type="V", mutation_class="S", vaf=None)
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True, include_somatic=0)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True, include_somatic=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "C")
        self.assertEqual(edit[3][3], "C")
        self.assertEqual(edit[3][4], "V")

        self.atoi_manual_transcript.edit('C', pos, mutation_type="V", mutation_class="S", vaf=None)
        edits, _ = self.atoi_manual_transcript.expressed_edits(include_rna_edits=True, include_somatic=0)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")
        edits, _ = self.atoi_manual_transcript.expressed_edits(include_rna_edits=True, include_somatic=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "C")
        self.assertEqual(edit[3][3], "C")
        self.assertEqual(edit[3][4], "V")

    def test_expressed_edit_with_overlapping_rna_germline_and_somatic_edits(self):
        pos = 9750162
        self.atoi_transcript.edit('C', pos, mutation_type="V", mutation_class="G", vaf=None)
        self.atoi_transcript.edit('A', pos, mutation_type="V", mutation_class="S", vaf=None)
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=2, include_somatic=0)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "C")
        self.assertEqual(edit[3][3], "C")
        self.assertEqual(edit[3][4], "V")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=0, include_somatic=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "RV")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=0, include_somatic=0)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=1, include_somatic=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "C")
        self.assertEqual(edit[3][2], "A")
        self.assertEqual(edit[3][3], "C")
        self.assertEqual(edit[3][4], "V")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=2, include_somatic=1)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][2], "C")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "RV")

        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=1, include_somatic=1)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")

        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True,
            include_germline=2, include_somatic=2)
        edit = edits[pos-1][0]
        self.assertEqual(edit[0], "I")
        self.assertEqual(edit[3][3], "I")
        self.assertEqual(edit[3][4], "R")

    def test_expressed_edit_with_overlapping_germline_and_somatic_edits(self):
        pos = 9750163
        self.atoi_transcript.edit('C', pos, mutation_type="V", mutation_class="G", vaf=None)
        self.atoi_transcript.edit('A', pos, mutation_type="V", mutation_class="S", vaf=None)
        edits, _ = self.atoi_transcript.expressed_edits(include_germline=1, include_somatic=1)
        self.assertNotIn(pos-1, edits)
        edits, _ = self.atoi_transcript.expressed_edits(include_germline=2, include_somatic=2)
        self.assertNotIn(pos-1, edits)

    def test_expressed_edit_with_deletion(self):
        pos = 9750162
        self.atoi_transcript.edit(1, pos, mutation_type="D")
        edits, _ = self.atoi_transcript.expressed_edits(include_rna_edits=True)
        self.assertNotIn(pos-1, edits)

    def test_seq_to_peptide_with_I_N(self):
        """checks whether sseq_to_peptide function can properly handle N and I"""
        seq = "AAIATIIGNIAA"
        pep, editing_positions, ambiguous_positions, peptide_warnings = transcript.seq_to_peptide(seq, return_positions=True)
        self.assertEqual(editing_positions, [0, 1.0, 2.0, 3.0])
        self.assertEqual(ambiguous_positions, [3.0])
        self.assertEqual(pep, "KMGE")

    def test_neopeptide_with_rna_edit(self):    
        self.transcript.edit("I", 5248244, mutation_type="R", mutation_class="R" )
        pep, protein = self.transcript.neopeptides(return_protein=True, include_rna_edits=True)
        self.assertEqual(protein[0:3], "MVR")


if __name__ == "__main__":
    unittest.main()
