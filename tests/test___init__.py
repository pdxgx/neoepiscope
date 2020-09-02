#!/usr/bin/env python
# coding=utf-8
"""
test_neoepiscope.py

Tests functions in neoepiscope.py.

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
from collections import defaultdict
import os.path as path, sys

from neoepiscope import *

import unittest
import filecmp
import os

neoepiscope_dir = os.path.dirname(
    os.path.dirname((os.path.abspath(getsourcefile(lambda: 0))))
)

def predicate(line):
    ''' whether reading first line of neoepiscope output '''
    if "Neoepiscope version" in line:
        return False
    return True

class TestGTFprocessing(unittest.TestCase):
    """Tests proper creation of dictionaries store GTF data"""

    def setUp(self):
        """Sets up gtf file and creates dictionaries for tests"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.gtf = os.path.join(self.base_dir, "Ychrom.gtf")
        self.gtf2 = os.path.join(self.base_dir, "Chr11.gtf")
        self.Ycds, self.Ytx = gtf_to_cds(self.gtf, "NA", pickle_it=False)
        self.Ytree = cds_to_tree(self.Ycds, "NA", pickle_it=False)
        self.cds11, self.tx11 = gtf_to_cds(self.gtf2, "NA", pickle_it=False)
        self.tree11 = cds_to_tree(self.cds11, "NA", pickle_it=False)
        self.lengths11 = cds_to_feature_length(self.cds11, self.tx11, "NA", pickle_it=False)
        self.counts = {'ENST00000325207.9_2': 1571.0, 'ENST00000325147.13_1': 372.0}
        self.tpm11 = feature_to_tpm_dict(self.counts, self.lengths11)

    def test_transcript_to_cds(self):
        """Fails if dictionary was built incorrectly"""
        self.assertEqual(len(self.Ycds.keys()), 220)
        start_test = [x for x in self.cds11['ENST00000429923.5_1'] if x[1] == "start_codon"]
        self.assertEqual(len(start_test), 1)
        self.assertEqual(start_test[0][2], 1891437)

    def test_cds_tree(self):
        """Fails if dictionary was built incorrectly"""
        self.assertEqual(len(self.Ytree.keys()), 1)
        self.assertEqual(len(self.Ytree["chrY"]), 2585)

    def test_transcript_extraction(self):
        """Fails if incorrect transcripts are pulled"""
        self.assertEqual(
            len(get_transcripts_from_tree("chrY", 150860, 150861, self.Ytree)), 10
        )
        self.coordinate_search = list(self.Ytree["chrY"].overlap(150860, 150861))
        self.transcripts = []
        for interval in self.coordinate_search:
            self.transcripts.append(interval[2])
        self.transcripts.sort()
        self.assertEqual(
            self.transcripts,
            [
                "ENST00000381657.7_3_PAR_Y",
                "ENST00000381663.8_3_PAR_Y",
                "ENST00000399012.6_3_PAR_Y",
                "ENST00000415337.6_3_PAR_Y",
                "ENST00000429181.6_2_PAR_Y",
                "ENST00000430923.7_3_PAR_Y",
                "ENST00000443019.6_2_PAR_Y",
                "ENST00000445062.6_2_PAR_Y",
                "ENST00000447472.6_3_PAR_Y",
                "ENST00000448477.6_2_PAR_Y",
            ],
        )

    def test_feature_lengths(self):
        """Fails if feature lengths are counted incorrectly"""
        self.assertEqual(self.lengths11['ENST00000332865.10_1'], 0.533)

    def test_tpm(self):
        """Fails if TPM is calculated incorrectly"""
        self.assertEqual(self.tpm11['ENST00000325207.9_2'], 820065.9484656778)


class TestVCFmerging(unittest.TestCase):
    """Tests proper merging of somatic and germline VCFS"""

    def setUp(self):
        """Sets up files to use for tests"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.varscan = os.path.join(self.base_dir, "Ychrom.varscan.vcf")
        self.germline = os.path.join(self.base_dir, "Ychrom.germline.vcf")
        self.precombined = os.path.join(self.base_dir, "Ychrom.combined.vcf")
        self.outvcf = os.path.join(self.base_dir, "Ychrom.testcombine.vcf")
        combine_vcf(self.germline, self.varscan, self.outvcf)

    def test_merge(self):
        """Fails if VCFs were merged improperly"""
        self.assertTrue(filecmp.cmp(self.outvcf, self.precombined))

    def tearDown(self):
        """Removes test file"""
        os.remove(self.outvcf)


class TestPrepHapCUT(unittest.TestCase):
    """Tests addition of unphased mutations to HapCUT2 output"""

    def setUp(self):
        """Sets up hapcut and vcf files to use for tests"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.hapcut = os.path.join(self.base_dir, "test.hapcut.out")
        self.vcf = os.path.join(self.base_dir, "test.vcf")
        self.phased_vcf = os.path.join(self.base_dir, "phased.vcf")
        self.germline_vcf = os.path.join(self.base_dir, "germline.vcf")
        self.complete_hapcut = os.path.join(self.base_dir, "complete_hapcut.out")
        self.test_hapcut = os.path.join(self.base_dir, "test_complete_hapcut.out")
        self.rbp_haplotypes = os.path.join(self.base_dir, "rbp.haplotypes")
        self.test_rbp = os.path.join(self.base_dir, "test.rbp.haplotypes")

    def test_haplotype_prep(self):
        """Tests that output of haplotype prep is correct for either regular
           hapcut output or phased VCFs
        """
        prep_hapcut_output(self.test_hapcut, self.hapcut, self.vcf)
        self.assertTrue(filecmp.cmp(self.test_hapcut, self.complete_hapcut))
        prep_hapcut_output(self.test_rbp, None, self.phased_vcf, phased_vcf=True)
        self.assertTrue(filecmp.cmp(self.rbp_haplotypes, self.test_rbp))


    def tearDown(self):
        """Removes test file"""
        os.remove(self.test_hapcut)
        os.remove(self.test_rbp)


class TestVAFpos(unittest.TestCase):
    """Tests fetching of VAF position from VCF file"""

    def setUp(self):
        """Sets up vcf files to use for tests"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.varscan = os.path.join(self.base_dir, "Ychrom.varscan.vcf")
        self.mutect = os.path.join(self.base_dir, "Ychrom.mutect.vcf")

    def test_position(self):
        """Fails if incorrect positions are returned"""
        self.assertEqual(get_vaf_pos(self.varscan), (5, 'FREQ'))
        self.assertEqual(get_vaf_pos(self.mutect), (4, 'FA'))


class TestHaplotypeProcessing(unittest.TestCase):
    """Tests proper processing of HAPCUT2 files"""

    def setUp(self):
        """Sets up input files and dictionaries to use for tests"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.ref_prefix = os.path.join(self.base_dir, "Chr11.ref")
        self.reference_index = bowtie_index.BowtieIndexReference(self.ref_prefix)
        self.Chr11gtf = os.path.join(self.base_dir, "Chr11.gtf")
        self.Chr11cds, self.Chr11tx = gtf_to_cds(self.Chr11gtf, "NA", pickle_it=False)
        for transcript in self.Chr11cds:
            for cds_block in self.Chr11cds[transcript]:
                cds_block[0] = cds_block[0].replace("chr", "")
        self.Chr11tree = cds_to_tree(self.Chr11cds, "NA", pickle_it=False)
        self.Chr11hapcut = os.path.join(self.base_dir, "Chr11.hapcut.out")
        self.rbp_ref_prefix = os.path.join(self.base_dir, "chr14_index")
        self.rbp_reference_index = bowtie_index.BowtieIndexReference(self.rbp_ref_prefix)
        self.Chr14gtf = os.path.join(self.base_dir, "Chr14.gtf")
        self.Chr14cds, self.Chr14tx = gtf_to_cds(self.Chr14gtf, "NA", pickle_it=False)
        for transcript in self.Chr14cds:
            for cds_block in self.Chr14cds[transcript]:
                cds_block[0] = cds_block[0].replace("chr", "")
        self.Chr14tree = cds_to_tree(self.Chr14cds, "NA", pickle_it=False)
        self.phased_hapcut = os.path.join(self.base_dir, "rbp.haplotypes")

    def test_hap_processing(self):
        """Fails if file is processed incorrectly"""
        Chr11_txs, homozygous_vars = process_haplotypes(self.Chr11hapcut, 
                                                        self.Chr11tree, 
                                                        phasing=True)
        phased_txs, phased_homozygous = process_haplotypes(self.phased_hapcut, 
                                                        self.Chr14tree,
                                                        phasing=True)
        self.assertEqual(sorted(Chr11_txs.keys()), ['ENST00000299106.8_2', 
                                                    'ENST00000398531.2_2', 
                                                    'ENST00000441717.3_2'])
        self.assertEqual(homozygous_vars["ENST00000299106.8_2"], 
                [
                    [
                        "11",
                        134018663,
                        "A",
                        "G",
                        "1",
                        "1",
                        "1/1:.:17:17:0:0%:17,0,0,0:.:2",
                        "V",
                    ]
                ]
        )
        self.assertEqual(Chr11_txs["ENST00000299106.8_2"],
            [
                [
                    [
                        "11",
                        134015873,
                        "GCAG",
                        4, 
                        "1",
                        "0",
                        "0/1:.:53:52:0:0%:22,30,0,0:.:2",
                        "D"
                    ],
                    [
                        "11",
                        134015876,
                        "",
                        "TT", 
                        "1",
                        "0",
                        "0/1:.:53:52:0:0%:22,30,0,0:.:2",
                        "I"
                    ]
                ],
                [
                    [
                        "11",
                        134019062,
                        "T",
                        "C",
                        "0",
                        "1",
                        "0/1:.:38:38:0:0%:31,7,0,0:.:2",
                        "V"
                    ]
                ]
            ]
        )
        self.assertEqual(
            Chr11_txs["ENST00000398531.2_2"],
            [
                [
                    [
                        "11",
                        71276862,
                        "GT",
                        2,
                        "0",
                        "1",
                        "0/1:.:53:52:0:0%:22,30,0,0:.:2",
                        "D",
                    ],
                    [
                        "11",
                        71276900,
                        "C",
                        "G",
                        "0",
                        "1",
                        "0/1:.:35:34:0:0%:19,15,0,0:.:2",
                        "V",
                    ],
                    [
                        "11",
                        71277000,
                        "",
                        "AA",
                        "0",
                        "1",
                        "0/1:.:35:34:0:0%:19,15,0,0:.:2",
                        "I",
                    ],
                ]
            ],
        )
        self.assertEqual(list(phased_txs.keys()), ['ENST00000409832.3'])
        self.assertEqual(phased_homozygous, {})
        self.assertEqual(
            phased_txs["ENST00000409832.3"],
            [
                [
                    [
                        "14",
                        19553372,
                        "G",
                        "A",
                        "0",
                        "1",
                        "0/1:647,136:783:99:19553372-1,19553372-2:1684,0,17385:GERMLINE*",
                        "V",
                    ],
                    [
                        "14",
                        19553436,
                        "C",
                        "T",
                        "1",
                        "0",
                        "0/1:740,103:843:99:19553372-2,19553372-1:1930,0,30239:3965.49:GERMLINE*",
                        "V",
                    ],
                    [
                        "14",
                        19553443,
                        "G",
                        "A",
                        "0",
                        "1",
                        "0/1:726,98:824:99:19553372-1,19553372-2:1889,0,29565:17731.95:GERMLINE*",
                        "V",
                    ],
                    [
                        "14",
                        19553764,
                        "A",
                        "G",
                        "0",
                        "1",
                        "0/1:1344,721:2065:99:19553372-1,19553372-2:15846,0,36781:726.04:GERMLINE*",
                        "V",
                    ],
                    [
                        "14",
                        19553795,
                        "G",
                        "A",
                        "0",
                        "1",
                        "0/1:16.22%:19553372-1,19553372-2:8761.31:SOMATIC",
                        "V",
                    ],
                ]
            ],
        )

    def test_maximum_clique(self):
        ht = [  ['11', 5246952, 'A', 'T', '0', '1', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'], 
                ['11', 5246956, 'G', 'A', '0', '1', 
                 '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'], 
                ['11', 5246956, 'G', 'T', '0', '1', 
                 '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'], 
                ['11', 5247812, 'A', 'T', '1', '0', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'],
                ['11', 5247832, 'AGCT', 4, '1', '0', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'D'],
                ['11', 5247834, 'CTT', 3, '1', '0', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'D'], 
                ['11', 5248161, '', 'A', '0', '1', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'I'],
                ['11', 5248161, '', 'T', '0', '1', 
                 '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'I'],                 
        ]
        cliques = list(transcript.get_haplotype_cliques(ht))
        for x in cliques:
            x.sort(key=itemgetter(1))
        sorted_cliques = sorted(cliques)
        self.assertEqual(
            sorted_cliques,
            [
                [
                 ('11', 5246952, 'A', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5246956, 'G', 'A', '0', '1', 
                  '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'),
                 ('11', 5248161, '', 'A', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'I')
                ],
                [
                 ('11', 5246952, 'A', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5246956, 'G', 'A', '0', '1', 
                  '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'),
                 ('11', 5248161, '', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'I')
                ],
                [
                 ('11', 5246952, 'A', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5246956, 'G', 'T', '0', '1', 
                  '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'),
                 ('11', 5248161, '', 'A', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'I')
                ],
                [
                 ('11', 5246952, 'A', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5246956, 'G', 'T', '0', '1', 
                  '0/1:.:53:52:0:3.0%:22,30,0,0:.:2', 'V'),
                 ('11', 5248161, '', 'T', '0', '1', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'I')
                ],
                [
                 ('11', 5247812, 'A', 'T', '1', '0', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5247832, 'AGCT', 4, '1', '0', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'D')
                ],
                [
                 ('11', 5247812, 'A', 'T', '1', '0', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2', 'V'),
                 ('11', 5247834, 'CTT', 3, '1', '0', 
                  '0/1:.:35:34:0:0.1%:19,15,0,0:.:2*', 'D')
                ],
            ]
        )

    def test_peptide_gathering(self):
        Chr11_txs = {
            "ENST00000398531.2_2": [
                [
                    [
                        "11",
                        71276651,
                        "CTC",
                        3,
                        "0",
                        "1",
                        "0/1:.:53:52:0:3.0%:22,30,0,0:.:2",
                        "D",
                    ],
                    [
                        "11",
                        71277229,
                        "A",
                        "C",
                        "0",
                        "1",
                        "0/1:.:35:34:0:15.7%:19,15,0,0:.:2",
                        "V",
                    ],
                    [
                        "11",
                        71277056,
                        "",
                        "AAA",
                        "0",
                        "1",
                        "0/1:.:35:34:0:0.1%:19,15,0,0:.:2",
                        "I",
                    ],
                ]
            ]
        }
        homozygous_vars = {'ENST00000299106.8_2': [
                                    [
                                        "11",
                                        134018663,
                                        "A",
                                        "G",
                                        "1",
                                        "1",
                                        "1/1:.:17:17:0:0%:17,0,0,0:.:2",
                                        "V",
                                    ]
                                ]
                        }
        transcript_blocks = self.Chr11cds["ENST00000398531.2_2"]
        neoepitopes, fasta = get_peptides_from_transcripts(
            Chr11_txs,
            homozygous_vars,
            (5, 'FREQ'),
            self.Chr11cds,
            True,
            False,
            False,
            self.reference_index,
            [8, 9, 10, 11],
            False,
            False,
            False,
            False,
            False,
            False,
            2,
            1,
            protein_fasta=True,
        )
        self.assertEqual(len(neoepitopes.keys()), 108)
        self.assertEqual(
            neoepitopes["CGCSQKCN"],
            [("11", 71277056, "", "AAA", "I", 0.001, "NA", "NA", "ENST00000398531.2_2")],
        )
        self.assertEqual(
            neoepitopes["PVCCPCKI"],
            [
                (
                    "11",
                    71277229,
                    "A",
                    "C",
                    "V",
                    0.157,
                    "PVCCQCKI",
                    "NA",
                    "ENST00000398531.2_2",
                )
            ],
        )
        self.assertEqual(
            neoepitopes["NKQDGESYE"], 
            [
                (
                    "11",
                    134018663,
                    "A",
                    "G",
                    "V",
                    0,
                    "NKQDGESYK",
                    "NA",
                    "ENST00000299106.8_2",
                )
            ]
        )
        self.assertEqual(sorted(neoepitopes.keys())[0], "CCGCGGCG")
        self.assertEqual(sorted(neoepitopes.keys())[-1], "YENPGKPDGVN")
        self.assertEqual(
            sorted(fasta["ENST00000398531.2_2"]),
            [
                "MGCCGCGGCGSGCGGCGSGCGGCGSGCGGYGSGCGGCGSSCCVPVCCCKPVCCCVPACSCSSCG"
                "SCGGSKGDCGSCGGSKGGCGSCGGSKGGCGSCGGSKGGCGSCGGSKGGCGSCGGSKGGCGS"
                "CGGSKGGCGSCGCSQKCNCCKPCCCSSGCGSCCQSSCCNPCCCQSSCCVPVCCQSSCCKPC"
                "CCQSSCCVPVCCPCKI"
            ],
        )

class TestExpression(unittest.TestCase):
    """Tests variant-level expression"""

    def setUp(self):
        """Sets up paths/variables"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.bam = os.path.join(self.base_dir, "test.rna.bam")
        self.neoepitopes = {
            'AAAAAAAAA': [
                (
                    '11', 
                    63401, 
                    'C', 
                    'T',
                    'V',
                    'NA', 
                    'AAAACAAAA', 
                    'NA', 
                    'NA', 
                    'TX1.2'
                )
            ]
        }
        self.ref_prefix = os.path.join(self.base_dir, "Chr11.ref")
        self.reference_index = bowtie_index.BowtieIndexReference(self.ref_prefix)

    def testSupport(self):
        """Test read support function"""
        expressed_vars, covered_vars = transcript_expression.get_expressed_variants(
                    self.bam, self.reference_index, self.neoepitopes
        )
        self.assertEqual(expressed_vars[('11', 63401, 'C', 'T', 'V')], 2)
        self.assertEqual(covered_vars[('11', 63401, 'C', 'T', 'V')], 4)
        self.assertEqual(len(expressed_vars.keys()), 1)
        self.assertEqual(len(covered_vars.keys()), 1)


class TestOutput(unittest.TestCase):
    """Tests function to write output"""

    def setUp(self):
        """Sets up paths and dictionaries"""
        self.base_dir = os.path.join(neoepiscope_dir, "tests")
        self.out_file = os.path.join(self.base_dir, "neoepiscope.out")
        self.correct_out = os.path.join(self.base_dir, "expected.neoepiscope.out")
        self.gtf = os.path.join(self.base_dir, "Chr11.gtf")
        self.cds, self.tx = gtf_to_cds(self.gtf, "NA", pickle_it=False)
        self.tools = {
            "netMHCpan4": ["netMHCpan", ["rank", "affinity"]],
            "netMHCIIpan3": ["netMHCIIpan", ["rank"]],
        }
        self.HLA_alleles = ["HLA*A01:01", "HLA*A02:01"]
        self.neoepitopes = {
            "CGCSQKCN": [
                (
                    "11",
                    71277056,
                    "",
                    "AAA",
                    "I",
                    0.001,
                    "NA",
                    "NA",
                    "ENST00000398531.2_2",
                    5,
                    10000.0,
                    1,
                    5,
                    150,
                    4,
                ),
                (
                    "11",
                    167789,
                    "A",
                    "T",
                    "V",
                    0.102,
                    "CGCSQCNN",
                    "NA",
                    "ENST00000410108.5_1",
                    0.5,
                    100,
                    0.5,
                    3,
                    150,
                    4,
                ),
            ],
            "PVCCPCKI": [
                (
                    "11",
                    71277229,
                    "A",
                    "C",
                    "V",
                    0.157,
                    "PVCCQCKI",
                    "NA",
                    "ENST00000398531.2_2",
                    10,
                    50.57,
                    1.2,
                    3,
                    10.1,
                    7,
                ),
                (
                    "11",
                    71277229,
                    "A",
                    "C",
                    "V",
                    0.203,
                    "PVCCQCKI",
                    "NA",
                    "ENST00000325113.8_1",
                    10,
                    50.57,
                    1.2,
                    3,
                    10.1,
                    7,
                ),
            ],
        }
        self.tpm_dict = {"ENST00000325113.8_1": 5.78, "ENST00000398531.2_2": 0.52}
        self.expressed_vars = defaultdict(int)
        self.expressed_vars[("11", 71277229, "A", "C", "V")] = 7
        self.expressed_vars[("11", 167789, "A", "T", "V")] = 5
        self.covered_vars = defaultdict(int)
        self.covered_vars[("11", 71277229, "A", "C", "V")] = 10
        self.covered_vars[("11", 167789, "A", "T", "V")] = 30

    def testwrite(self):
        """Tests that output file is written correctly"""

        from sys import version_info
        write_results(self.out_file, self.HLA_alleles, self.neoepitopes, 
                      self.tools, self.tx, self.tpm_dict, None, self.expressed_vars,
                      self.covered_vars)
        if version_info[0] < 3:
            from itertools import izip, ifilter
            with open(self.out_file) as fh1, open(self.correct_out) as fh2:
                f1 = ifilter(predicate, fh1)
                f2 = ifilter(predicate, fh2)
                test = all(x == y for x, y in izip(f1, f2))
        else:
            with open(self.out_file) as fh1, open(self.correct_out) as fh2:
                f1 = filter(predicate, fh1)
                f2 = filter(predicate, fh2)
                test = all(x == y for x, y in zip(f1, f2))
        self.assertTrue(test)

    def tearDown(self):
        """Removes test file"""
        os.remove(self.out_file)


if __name__ == "__main__":
    unittest.main()
