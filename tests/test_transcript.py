import unittest
from neoepiscope.transcript import *
unittest.TestCase.maxDiff = None
import os

class TestTranscript(unittest.TestCase):
    """Tests transcript object construction"""
    def setUp(self):
        """Sets up gtf file and creates dictionaries for tests"""
        self.gtf = os.path.join(
                        os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 
                        'tests', 'Chr11.gtf'
                        )
        self.cds = gtf_to_cds(self.gtf, 'NA', pickle_it=False)
        self.ref_prefix = os.path.join(
                       os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 
                       'tests', 'Chr11.ref'
                    )
        self.reference_index = bowtie_index.BowtieIndexReference(
                                                    self.ref_prefix)
        self.transcript = Transcript(self.reference_index, 
                                    [[str(chrom).replace('chr', ''), 'blah', seq_type,
                                      str(start), str(end), '.', 
                                      strand] for (chrom, seq_type, start, 
                                                    end, strand) in 
                                      self.cds['ENST00000335295.4_1']],
                                      'ENST00000335295.4_1')
        self.fwd_transcript = Transcript(self.reference_index, 
                                    [[str(chrom).replace('chr', ''), 'blah', seq_type,
                                      str(start), str(end), '.', 
                                      strand] for (chrom, seq_type, start, 
                                                    end, strand) in 
                                      self.cds['ENST00000308020.5_1']],
                                      'ENST00000308020.5_1')
        self.all_coding_transcript = Transcript(self.reference_index, 
                                    [[str(chrom).replace('chr', ''), 'blah', seq_type,
                                      str(start), str(end), '.', 
                                      strand] for (chrom, seq_type, start, 
                                                    end, strand) in 
                                      self.cds['ENST00000317078.1_1']],
                                      'ENST00000317078.1_1')
    def test_transcript_structure(self):
        """Fails if structure of unedited transcript is incorrect"""
        self.assertEqual(len(self.transcript.annotated_seq()), 3)
        self.assertEqual(len(self.transcript.annotated_seq()[0][0]), 142)
        self.assertEqual(len(self.transcript.annotated_seq()[1][0]), 223)
        self.assertEqual(len(self.transcript.annotated_seq()[2][0]), 263)
        self.assertEqual(self.transcript.annotated_seq()[0][1], 'R')
        self.assertEqual(self.transcript.intervals, [5246692, 5246955, 
                                                     5247805, 5248028,
                                                     5248158, 5248300])
        self.assertEqual(self.transcript.start_codon, 5248249)
        self.assertEqual(self.transcript._start_codon, 5248248)
        self.assertEqual(self.transcript.stop_codon, 5246828)
        self.assertEqual(self.transcript._stop_codon, 5246827)
        self.assertTrue(self.transcript.rev_strand)
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(self.transcript.deletion_intervals, [])
        self.assertEqual(self.transcript.transcript_id, 'ENST00000335295.4_1')
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
    # Transcript sequence editing tests
    def test_irrelevant_edit(self):
        """Fails if edit is made for non-exon position"""
        self.transcript.edit('A', 5248155)
        relevant_edits = self.transcript.expressed_edits()
        self.assertEqual(self.transcript.edits[5248154], [('A', 'V', 'S', 
                                                          ('11', 5248155, 
                                                           'C', 'A',
                                                           'V', None))])
        self.assertEqual(relevant_edits[0], {})
        self.assertEqual(relevant_edits[1], [(5246692, 'R', ()),
                                             (5246955, 'R', ()),
                                             (5247805, 'R', ()), 
                                             (5248028, 'R', ()), 
                                             (5248158, 'R', ()), 
                                             (5248300, 'R', ())])
        self.assertEqual(len(self.transcript.annotated_seq()), 3)
        self.assertEqual(len([x for x in self.transcript.annotated_seq() 
                                        if x[1] != 'R']), 0)
    def test_relevant_edit(self):
        """Fails if edit is not made for position within exon"""
        self.transcript.edit('A', 5248299)
        relevant_edits = self.transcript.expressed_edits()
        self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 'S',
                                                          ('11', 5248299, 
                                                            'T', 'A', 
                                                            'V', None))])
        self.assertEqual(relevant_edits[0][5248298], [('A', 'V', 'S',
                                                      ('11', 5248299, 'T', 
                                                        'A', 'V', None))])
        self.assertEqual(relevant_edits[1], [(5246692, 'R', ()),
                                             (5246955, 'R', ()),
                                             (5247805, 'R', ()), 
                                             (5248028, 'R', ()), 
                                             (5248158, 'R', ()), 
                                             (5248300, 'R', ())])
        self.assertEqual(len(self.transcript.annotated_seq()), 5)
        self.assertEqual(self.transcript.annotated_seq()[1][0], 'T')
        self.assertEqual(self.transcript.annotated_seq()[1][1], 'S')
        self.assertEqual(len([x for x in self.transcript.annotated_seq() 
                                        if x[1] != 'R']), 1)
    def test_reset_to_reference(self):
        """Fails if transcript is not reset to reference"""
        self.transcript.edit('A', 5248299)
        self.transcript.reset(reference=True)
        self.assertEqual(self.transcript.edits, {})
    def test_edit_and_save(self):
        """Fails if edits aren't saved"""
        self.transcript.edit('A', 5248299)
        self.transcript.edit(3, 5246694, mutation_type='D')
        self.transcript.save()
        self.assertEqual(self.transcript.last_edits[5248298], [('A', 'V', 
                                                                'S', 
                                                               ('11', 
                                                                5248299, 
                                                                'T', 'A', 
                                                                'V', 
                                                                None))])
        self.assertEqual(self.transcript.last_deletion_intervals,
                            [(5246692, 5246695, 'S', ('11', 5246694, 'TTG',
                                                        '', 'D', None))])
    def test_reset_to_save_point(self):
        """Fails if new edit not erased or old edits not retained"""
        self.transcript.edit('A', 5248299)
        self.transcript.edit(3, 5246694, mutation_type='D')
        self.transcript.save()
        self.transcript.edit('G', 5248165)
        self.assertEqual(self.transcript.edits[5248164], [('G', 'V', 'S',
                                                          ('11', 5248165, 
                                                            'C', 'G',
                                                            'V', None))])
        self.transcript.reset(reference=False)
        self.assertNotIn(2182387, self.transcript.edits)
        self.assertEqual(self.transcript.last_edits[5248298], [('A', 'V', 
                                                                'S', 
                                                               ('11',
                                                                5248299, 
                                                                'T', 'A', 
                                                                'V', 
                                                                None))])
        self.assertEqual(self.transcript.last_deletion_intervals,
                            [(5246692, 5246695, 'S', ('11', 5246694, 'TTG', 
                                                       '', 'D', None))])
        self.assertNotEqual(self.transcript.edits, {})
    def test_SNV_seq(self):
        """Fails if SNV is edited incorrectly"""
        self.transcript.edit('A', 5248299)
        self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 
                                                                'S', 
                                                               ('11', 
                                                                5248299, 
                                                                'T', 'A',
                                                                'V',
                                                                None))])
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(seq[0], ('AC', 'R', [()], 5248301))
        self.assertEqual(seq[1], ('T', 'S', [('11', 5248299, 'T', 'A', 
                                                'V', None)], 
                                    5248299))
        self.assertEqual(len(seq[2][0]), 139)
        self.assertEqual(len(seq[3][0]), 223)
        self.assertEqual(len(seq[4][0]), 263)
    def test_inside_insertion(self):
        """Fails if indel within exon is inserted incorrectly"""
        self.transcript.edit('Q', 5248165, mutation_type='I')
        self.assertEqual(self.transcript.edits[5248164], [('Q', 'I', 'S',
                                                            ('11', 5248165, 
                                                             '', 'Q',
                                                             'I', None))])
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 6)
        self.assertEqual(seq[1], ('Q', 'S', [('11', 5248165, '', 
                                              'Q', 'I', None)], 
                                  5248165))
        self.assertEqual(len(seq[0][0]), 136)
        self.assertEqual(len(seq[2][0]), 1)
        self.assertEqual(len(seq[3][0]), 5)
        self.assertEqual(len(seq[4][0]), 223)
        self.assertEqual(len(seq[5][0]), 263)
    def test_adjacent_indel(self):
        """Fails if indel right before exon is inserted incorrectly"""
        self.transcript.edit('Q', 5248029, mutation_type='I')
        self.assertEqual(self.transcript.edits[5248028], [('Q', 'I', 'S',
                                                            ('11', 5248029, 
                                                                '', 'Q',
                                                            'I', None))])
        self.assertEqual(self.transcript.deletion_intervals, [])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(seq[1], ('Q', 'S', [('11', 5248029, '', 'Q', 'I', 
                                              None)], 
                                    5248029)) 
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[2][0]), 1)
        self.assertEqual(len(seq[3][0]), 222)
        self.assertEqual(len(seq[4][0]), 263)
    def test_inside_deletion(self):
        """Fails if deletion completely within an exon is made improperly"""
        self.transcript.edit(3, 5246700, mutation_type='D')
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(self.transcript.deletion_intervals, [(5246698, 
                                                           5246701, 'S',
                                                           ('11', 
                                                            5246700, 
                                                            'TGA', '',
                                                            'D',
                                                            None))])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 5)
        self.assertEqual(seq[3], ('', 'S', [('11', 5246700, 'TGA', '', 'D', 
                                            None)],
                                    5246700))
        self.assertEqual(seq[4], ('TTGCAA', 'R', [()], 5246699))
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[1][0]), 223)
        self.assertEqual(len(seq[2][0]), 254)
    def test_overlapping_deletion(self):
        """Fails if deletion overlapping a junction is incorrect"""
        self.transcript.edit(10, 5246950, mutation_type='D')
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(self.transcript.deletion_intervals, [(5246948,
                                                           5246958, 'S',
                                                           ('11', 
                                                            5246950,
                                                            'CCAGGAGCTG',
                                                            '', 'D', 
                                                            None))])
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 4)
        self.assertEqual(seq[2], ('', 'S', [('11', 5246950,'CCAGGAGCTG', 
                                             '', 'D', None)], 5246950))
        self.assertEqual(len(seq[0][0]), 142)
        self.assertEqual(len(seq[1][0]), 223)
        self.assertEqual(len(seq[3][0]), 256)
    def test_spanning_deletion(self):
        """Fails if deletion spanning two exons is incorrect"""
        self.transcript.edit(137, 5248025, mutation_type='D')
        self.assertEqual(self.transcript.edits, {})
        self.assertEqual(self.transcript.deletion_intervals[0][0:3],
                                                (5248023, 5248160, 'S'))
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 4)
        self.assertEqual(len(seq[1][2][0][2]), 137)
        self.assertEqual(seq[1][0:2], ('', 'S'))
        self.assertEqual(len(seq[0][0]), 140)
        self.assertEqual(len(seq[2][0]), 218)
        self.assertEqual(len(seq[3][0]), 263)
    def test_deletion_over_transcript_start(self):
        """Fails if deletion spanning start of transcript is incorrect"""
        self.transcript.edit(5, 5246692, mutation_type='D')
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 4)
        self.assertEqual(seq[-1][0], '')
        self.assertEqual(len(seq[-2][0]), 260)
        self.fwd_transcript.edit(10, 450275, mutation_type='D')
        seq2 = self.fwd_transcript.annotated_seq()
        self.assertEqual(len(seq2), 13)
        self.assertEqual(seq2[0][0], '')
        self.assertEqual(len(seq2[1][0]), 353)
    def test_deletion_of_transcript(self):
        """Fails if deletion of entire transcript is incorrect"""
        self.transcript.edit(1700, 5246690, mutation_type='D')
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 1)
        self.assertEqual(seq[0][0], '')
        self.assertEqual(len(seq[0][2][0][2]), 1700)
    def test_compound_variants(self):
        """Fails if transcript with multiple variant types is incorrect"""
        self.transcript.edit(137, 5248025, mutation_type='D')
        self.transcript.edit('Q', 5248165, mutation_type='I')
        self.transcript.edit('A', 5248299)
        self.assertEqual(len(self.transcript.edits.keys()), 2)
        self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 'S', 
                                                           ('11', 5248299, 
                                                            'T', 'A','V', 
                                                                None))]) 
        self.assertEqual(self.transcript.edits[5248164], [('Q', 'I', 'S',
                                                            ('11', 5248165, 
                                                             '', 'Q',
                                                            'I', None))])
        self.assertEqual(self.transcript.deletion_intervals[0][0:3],
                                                (5248023, 5248160, 'S'))
        seq = self.transcript.annotated_seq()
        self.assertEqual(len(seq), 9)
        self.assertEqual(seq[0], ('AC', 'R', [()], 5248301))
        self.assertEqual(seq[1], ('T', 'S', [('11', 5248299, 'T', 'A', 'V', 
                                                None)], 
                                    5248299))
        self.assertEqual(len(seq[2][0]), 133)
        self.assertEqual(seq[3], ('Q', 'S', [('11', 5248165, '', 'Q', 'I', 
                                                None)], 
                                    5248165))
        self.assertEqual(seq[4], ('G', 'R', [()], 5248165))
        self.assertEqual(seq[5], ('GGC', 'R', [()], 5248164))
        self.assertEqual(len(seq[6][2][0][2]), 137)
        self.assertEqual(len(seq[7][0]), 218)
        self.assertEqual(len(seq[8][0]), 263)
    # Neopeptide tests
    def test_no_mutations_peptides(self):
        """Fails if peptides are returned for unmutated sequence"""
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(peptides, [])
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(rev_peptides, [])
    def test_noncoding_mutation_peptides(self):
        """Fails if peptides are returned for mutation in noncoding 
            sequence"""
        self.fwd_transcript.edit('G', 450286)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(peptides, [])
        self.transcript.edit('A', 5248266)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(rev_peptides, [])
    def test_synonymous_snv_peptides(self):
        """Fails if peptides are returned for a synonymous snv"""
        self.fwd_transcript.edit('A', 450464)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(peptides, [])
        self.transcript.edit('A', 5248005)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(rev_peptides, [])
    def test_missense_snv_peptides(self):
        """Fails if incorrect peptides are returned for missense SNV"""
        self.fwd_transcript.edit('T', 450502)
        peptides = self.fwd_transcript.neopeptides().keys()
        F_peptides = [pep for pep in peptides if 'F' in pep]
        self.assertEqual(len(peptides), 38)
        self.assertEqual(len(peptides), len(F_peptides))
        self.assertEqual(sorted(peptides)[0], 'AGGPRPEF')
        self.assertEqual(sorted(peptides)[-1], 'RRDAGGPRPEF')
        self.transcript.edit('T', 5248006)
        rev_peptides = self.transcript.neopeptides().keys()
        N_peptides = [pep for pep in rev_peptides if 'N' in pep]
        self.assertEqual(len(rev_peptides), 38)
        self.assertEqual(len(rev_peptides), len(N_peptides))
        self.assertEqual(sorted(rev_peptides)[0], 'GRLLVVYPWN')
        self.assertEqual(sorted(rev_peptides)[-1], 'YPWNQRFFESF')
    def test_in_frame_insertion_peptides(self):
        """Fails if incorrect peptides are returned for in-frame 
            insertion"""
        self.fwd_transcript.edit('AAA', 450551, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides().keys()
        K_peptides = [pep for pep in peptides if 'K' in pep]
        self.assertEqual(len(peptides), 38)
        self.assertEqual(len(peptides), len(K_peptides))
        self.assertEqual(sorted(peptides)[0], 'ASLEEPPDGPK')
        self.assertEqual(sorted(peptides)[-1], 'SLEEPPDGPKS')
        self.transcript.edit('TTT', 5247986, mutation_type='I')
        rev_peptides = self.transcript.neopeptides().keys()
        K_peptides = [pep for pep in peptides if 'K' in pep]
        self.assertEqual(len(rev_peptides), 38)
        self.assertEqual(len(rev_peptides), len(K_peptides))
        self.assertEqual(sorted(rev_peptides)[0], 'ESKFGDLS')
        self.assertEqual(sorted(rev_peptides)[-1], 'YPWTQRFFESK')
    def test_synonymous_inframe_insertion_peptides(self):
        """Fails if incorrect peptides are returned for in insertion into
            a codon that maintains the AA sequence of that codon"""
        self.fwd_transcript.edit('AAA', 450502, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 38)
        self.transcript.edit('TTT', 5246874, mutation_type='I')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 34)
    def test_frameshift_insertion(self):
        """Fails if incorrect peptides are returned for frameshift
            insertion"""
        self.fwd_transcript.edit('AAAAA', 473925, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 108)
        self.assertEqual(sorted(peptides)[0], 'ASILVFLK')
        self.assertEqual(sorted(peptides)[-1], 'VLESHKLKTGH')
        self.transcript.edit('AAAAA', 5247883, mutation_type='I')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(sorted(rev_peptides)[0], 'AFSDGLAHLV')
        self.assertEqual(sorted(rev_peptides)[-1], 'VFTTSRAPLPH')
    def test_in_frame_deletion_peptides(self):
        """Fails if incorrect peptides are given for in-frame deletion"""
        self.fwd_transcript.edit(3, 450555, mutation_type='D')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 34)
        self.assertEqual(sorted(peptides)[0], 'DGPSGQAT')
        self.assertEqual(sorted(peptides)[-1], 'SLEEPPDGPSG')
        self.transcript.edit(3, 5247858, mutation_type='D')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 34)
        self.assertEqual(sorted(rev_peptides)[0], 'ALSELHCD')
        self.assertEqual(sorted(rev_peptides)[-1], 'TFALSELHCDK')
    def test_synonymous_inframe_deletion_peptides(self):
        """Fails if incorrect peptides are returned for in insertion into
            a codon that maintains the AA sequence of that codon"""
        self.fwd_transcript.edit(3, 473918, mutation_type='D')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 34)
        self.transcript.edit(3, 5247922, mutation_type='D')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 30)
    def test_frameshift_deletion(self):
        """Fails if incorrect peptides are returned for frameshift
            deletion"""
        self.fwd_transcript.edit(5, 473912, mutation_type='D')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 40)
        self.assertEqual(sorted(peptides)[0], 'ASSFLMFW')
        self.assertEqual(sorted(peptides)[-1], 'YNTKRGIVASS')
        self.transcript.edit(5, 5247930, mutation_type='D')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 32)
        self.assertEqual(sorted(rev_peptides)[0], 'AVMGNPKVKG')
        self.assertEqual(sorted(rev_peptides)[-1], 'VMGNPKVKGQE')
    def test_nonstop_mutation_peptides(self):
        """Fails if mutation altering stop codon does not return peptides
            past the end of the original peptide to the new stop"""
        self.fwd_transcript.edit('A', 490580)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 768)
        self.assertEqual(sorted(peptides)[0], 'AARVGARP')
        self.assertEqual(sorted(peptides)[-1], 'YTGSTSTSPAA')
        self.transcript.edit('G', 5246828)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 84)
        self.assertEqual(sorted(rev_peptides)[0], 'AHKYHYAR')
        self.assertEqual(sorted(rev_peptides)[-1], 'YHYARFLAVQF')
    def test_split_start(self):
        """Fails if split start codon is handled improperly"""
        self.fwd_transcript.edit('AT', 450419, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 278)
        self.assertEqual(sorted(peptides)[0], 'AAAAPSPR')
        self.assertEqual(sorted(peptides)[-1], 'WRSRLTGRLPA')
        self.transcript.edit('T', 5248290, mutation_type='I')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 42)
        self.assertEqual(sorted(rev_peptides)[0], 'ATSNRHHG')
        self.assertEqual(sorted(rev_peptides)[-1], 'TSNRHHGASDS')
    def test_start_lost_peptides(self):
        """Fails if mutation altering start codon does not return peptides
           from a new start codon"""
        self.fwd_transcript.edit('T', 450456)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(peptides, [])
        # Next start immediately followed by stop for fwd strand transcript
        self.transcript.edit('G', 5248251)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 122)
        self.assertEqual(sorted(rev_peptides)[0], 'AGCWWSTL')
        self.assertEqual(sorted(rev_peptides)[-1], 'WWSTLGPRGSL')
    def test_start_lost_and_new_inframe_start(self):
        """Fails if peptides aren't returned from a new in frame start 
           codon when the original is disrupted"""
        self.fwd_transcript.edit('ATG', 450446, mutation_type='I')
        self.fwd_transcript.edit('T', 450456)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 20)
        ### Got to fix this one below
        self.transcript.edit('CAT', 5248263, mutation_type='I')
        self.transcript.edit('G', 5248251)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 24)
    def test_start_lost_and_new_out_of_frame_start(self):
        """Fails if peptides aren't returned from a new out of frame start 
            codon when the original is disrupted"""
        self.fwd_transcript.edit('ATG', 450445, mutation_type='I')
        self.fwd_transcript.edit('T', 450456)
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 98)
        self.transcript.edit('CAT', 5248280, mutation_type='I')
        self.transcript.edit('G', 5248251)
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 22)
    def test_skipping_new_start(self):
        """Fails if peptides are returned from a new start codon when the 
            original is retained"""
        self.fwd_transcript.edit('ATG', 450445, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides(
                            only_downstream=True
                            ).keys()
        self.assertEqual(peptides, [])
        self.transcript.edit('CAT', 5248280, mutation_type='I')
        rev_peptides = self.transcript.neopeptides(
                            only_downstream=True
                            ).keys()
        self.assertEqual(rev_peptides, [])
    def test_compound_indel_peptides(self):
        """Fails if incorrect peptides are returned with complementary
            indels are introduced (i.e. frameshift then return to frame)"""
        self.fwd_transcript.edit(4, 473924, mutation_type='D')
        self.fwd_transcript.edit('AAAA', 473952, mutation_type='I')
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 74)
        self.assertEqual(sorted(peptides)[0], 'ASILVFFL')
        self.assertEqual(sorted(peptides)[-1], 'VFFLESHKLKT')
        self.transcript.edit(4, 5247921, mutation_type='D')
        self.transcript.edit('AAAA', 5247933, mutation_type='I')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 50)
    def test_all_coding_tx(self):
        """Fails if transcript that is all coding sequence is 
            handled improperly"""
        self.all_coding_transcript.edit('C', 5810036)
        peptides = self.all_coding_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 16)
    def test_germline_vs_somatic(self):
        """"Fails if incorrect peptides are returned for different
            germline/somatic mutation inclusion decisions"""
        self.fwd_transcript.edit('T', 450502)
        self.fwd_transcript.edit('T', 450503, mutation_class='G')
        self.assertEqual(len(self.fwd_transcript.neopeptides().keys()), 38)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_germline=1).keys()), 38)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=0).keys()), 0)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=1).keys()), 0)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=2).keys()), 0)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_somatic=2,
                                        include_germline=0).keys()), 0)
        self.assertEqual(len(self.fwd_transcript.neopeptides(
                                        include_somatic=2,
                                        include_germline=1).keys()), 0)
        self.transcript.edit('T', 5248006)
        self.transcript.edit('C', 5248007, mutation_class='G')
        self.assertEqual(len(self.transcript.neopeptides().keys()), 38)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_germline=1).keys()), 38)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=0).keys()), 0)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=1).keys()), 38)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_somatic=0,
                                        include_germline=2).keys()), 0)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_somatic=2,
                                        include_germline=0).keys()), 0)
        self.assertEqual(len(self.transcript.neopeptides(
                                        include_somatic=2,
                                        include_germline=1).keys()), 38)
    def test_compound_all(self):
        """Fails if incorrect peptides are returned when multiple
            germline/somatic mutations are introduced"""
        self.fwd_transcript.edit('AAA', 490579, 
            mutation_type='I', mutation_class='G')
        self.fwd_transcript.edit('C', 490580, mutation_type='V')
        self.fwd_transcript.edit('A', 490582, mutation_type='I')
        self.fwd_transcript.edit('A', 490586, mutation_type='V')
        self.fwd_transcript.edit('G', 490587, mutation_type='D')
        self.fwd_transcript.edit('C', 490588, mutation_type='V')
        self.fwd_transcript.edit('A', 490593, 
            mutation_type='V', mutation_class='G')  
        self.fwd_transcript.edit("GAGGAGGAGGAG",450536,mutation_type="I")
        self.fwd_transcript.save()
        self.fwd_transcript.edit('T', 450457, mutation_type='D')
        self.fwd_transcript.reset()                                 
        peptides = self.fwd_transcript.neopeptides().keys()
        self.assertEqual(len(peptides), 58)
        self.assertEqual(sorted(peptides)[5], 'APTPNKRTYP')
        self.assertEqual(sorted(peptides)[-1], 'VPAGRASLEEE')
        self.assertEqual(sorted(peptides)[-5], 'SLEEEEEEP')
        peptides = self.fwd_transcript.neopeptides(
            include_germline=1).keys()
        self.assertEqual(len(peptides), 62)
        self.assertEqual(sorted(peptides)[0], 'AEGEGAPTPNK')
        self.assertEqual(sorted(peptides)[32], 'EGEGAPTPNKR')
        self.assertEqual(sorted(peptides)[-5], 'SLEEEEEEP')
        # Reverse transcript      
        self.transcript.edit(4, 5247921, mutation_type='D')
        self.transcript.edit('AAAA', 5247933, mutation_type='I')
        self.transcript.edit(3, 5248170, mutation_type='D', 
            mutation_class='G')
        self.transcript.edit('T', 5246872, mutation_type='V')
        rev_peptides = self.transcript.neopeptides().keys()
        self.assertEqual(len(rev_peptides), 88)
        self.assertEqual(sorted(rev_peptides)[0], 'AAYQKMVA')
        self.assertEqual(sorted(rev_peptides)[-1], 'YQKMVAGVANA')
        rev_peptides = self.transcript.neopeptides(
                        include_germline=1).keys()
        self.assertEqual(len(rev_peptides), 122)            
        self.assertEqual(sorted(rev_peptides)[13], 'DEVGGALG')
        self.assertEqual(sorted(rev_peptides)[-13], 'VNVDEVGGALG')
        self.assertEqual(sorted(rev_peptides)[0], 'AAYQKMVA')

if __name__ == '__main__':
    unittest.main()