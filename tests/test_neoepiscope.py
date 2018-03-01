from neoepiscope import *

import unittest
import filecmp
import os 
neoepiscope_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

class TestGTFprocessing(unittest.TestCase):
	"""Tests proper creation of dictionaries store GTF data"""
	def setUp(self):
		"""Sets up gtf file and creates dictionaries for tests"""
		self.base_dir = neoepiscope_dir
		self.gtf = os.path.join(self.base_dir
							, 'tests', 'Ychrom.gtf'
						)
		self.Ycds = gtf_to_cds(self.gtf, 'NA', pickle_it=False)
		self.Ytree = cds_to_tree(self.Ycds, 'NA', pickle_it=False)
	def test_transcript_to_CDS(self):
		"""Fails if dictionary was built incorrectly"""
		self.assertEqual(len(self.Ycds.keys()), 164)
	def test_CDS_tree(self):
		"""Fails if dictionary was built incorrectly"""
		self.assertEqual(len(self.Ytree.keys()), 1)
		self.assertEqual(len(self.Ytree['Y']), 2174)
	def test_transcript_extraction(self):
		"""Fails if incorrect transcripts are pulled"""
		self.assertEqual(len(get_transcripts_from_tree(
														'Y', 150860, 
														150861,
														self.Ytree)), 
														3
														)
		self.coordinate_search = list(self.Ytree['Y'].search(150860,
															 150861))
		self.transcripts = []
		for interval in self.coordinate_search:
			self.transcripts.append(interval[2])
		self.transcripts.sort()
		self.assertEqual(
						 self.transcripts, ['ENST00000381657.7_3_PAR_Y',
											'ENST00000381663.8_3_PAR_Y',
											'ENST00000399012.6_3_PAR_Y']
						)
class TestVCFmerging(unittest.TestCase):
	"""Tests proper merging of somatic and germline VCFS"""
	def setUp(self):
		"""Sets up files to use for tests"""
		self.base_dir = neoepiscope_dir
		self.varscan = os.path.join(self.base_dir, 'tests', 'Ychrom.varscan.vcf'
						)                
		self.germline = os.path.join(self.base_dir, 'tests', 'Ychrom.germline.vcf'
						)   
		self.precombined = os.path.join(self.base_dir, 'tests', 'Ychrom.combined.vcf'
						)  
		self.outvcf = os.path.join(self.base_dir, 'tests', 'Ychrom.testcombine.vcf'
						)  
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
		self.base_dir = neoepiscope_dir
		self.hapcut = os.path.join(self.base_dir, 'tests', 'test.hapcut.out'
						)
		self.vcf = os.path.join(self.base_dir, 'tests', 'test.vcf'
						)
		self.complete_hapcut = os.path.join(self.base_dir, 'tests', 'complete_hapcut.out'
						)
		self.test_hapcut = os.path.join(self.base_dir, 'tests', 'test_complete_hapcut.out'
						)
	def test_haplotype_prep(self):
		"""Tests that output of haplotype prep is correct"""
		prep_hapcut_output(self.test_hapcut, self.hapcut, self.vcf)
		self.assertTrue(filecmp.cmp(self.test_hapcut,
									self.complete_hapcut))
	def tearDown(self):
		"""Removes test file"""
		os.remove(self.test_hapcut)
class TestVAFpos(unittest.TestCase):
	"""Tests fetching of VAF position from VCF file"""
	def setUp(self):
		"""Sets up vcf files to use for tests"""
		self.base_dir = neoepiscope_dir
		self.varscan = os.path.join(self.base_dir, 'tests', 'Ychrom.varscan.vcf'
						)
		self.mutect = os.path.join(self.base_dir, 'tests', 'Ychrom.mutect.vcf'
						)
	def test_position(self):
		"""Fails if incorrect positions are returned"""
		self.assertEqual(get_VAF_pos(self.varscan), 5)
		self.assertEqual(get_VAF_pos(self.mutect), None)
class TestHaplotypeProcessing(unittest.TestCase):
	"""Tests proper processing of HAPCUT2 files"""
	def setUp(self):
		"""Sets up input files and dictionaries to use for tests"""
		self.base_dir = neoepiscope_dir
		self.ref_prefix = os.path.join(
					self.base_dir, 'tests', 'Chr11.ref'
				)
		self.reference_index = bowtie_index.BowtieIndexReference(
														self.ref_prefix)
		self.Chr11gtf = os.path.join(self.base_dir, 'tests', 'Chr11.gtf'
						)
		self.Chr11cds = gtf_to_cds(self.Chr11gtf, 'NA', pickle_it=False)
		self.Chr11tree = cds_to_tree(self.Chr11cds, 'NA',
										pickle_it=False)
		self.Chr11hapcut = os.path.join(self.base_dir, 'tests', 'Chr11.hapcut.out'
						)
	def test_hap_processing(self):
		"""Fails if file is processed incorrectly"""
		Chr11_txs = process_haplotypes(self.Chr11hapcut, self.Chr11tree)
		self.assertEqual(sorted(Chr11_txs.keys()), 
								['ENST00000398531.2_2'])
		self.assertEqual(Chr11_txs['ENST00000398531.2_2'],
						[[['11', 71276862, 'TGT', 2, '0', '0', 
						   '0/0:.:53:52:0:0%:22,30,0,0:.:2', 'D'], 
						  ['11', 71276900, 'C', 'G', '0', '0',
						   '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'V'],
						  ['11', 71277000, 'G', 'AA', '0', '0',
						   '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'I']]]
			)
	def test_peptide_gathering(self):
		Chr11_txs = {'ENST00000398531.2_2': [
							[['11', 71276651, 'CTC', 3, '0', '1', 
							  '0/0:.:53:52:0:3.0%:22,30,0,0:.:2', 'D'], 
							 ['11', 71277229, 'A', 'C', '0', '1',
							  '0/0:.:35:34:0:15.7%:19,15,0,0:.:2', 'V'],
							 ['11', 71277056, 'G', 'AAA', '1', '1',
							  '0/0:.:35:34:0:0.1%:19,15,0,0:.:2', 'I']]
							 ]
					}
		transcript_blocks = self.Chr11cds['ENST00000398531.2_2']
		neoepitopes = get_peptides_from_transcripts(Chr11_txs, 5, 
													self.Chr11cds,
													True, False, False, 
													self.reference_index,
													[8,9,10,11])
		self.assertEqual(len(neoepitopes.keys()), 70)
		self.assertEqual(neoepitopes['CGCSQKCN'], [('11', 71277056, '',
												'AAA', 'I', 0.1, 
												'ENST00000398531.2_2')])
		self.assertEqual(neoepitopes['PVCCPCKI'], [('11', 71277229, 
												'A', 'C', 'V', 15.7, 
												'ENST00000398531.2_2')])
		self.assertEqual(sorted(neoepitopes.keys())[0], 'CCGCGGCG')
		self.assertEqual(sorted(neoepitopes.keys())[-1], 'VPVCCPCKI')
class TestBindingPrediction(unittest.TestCase):
	"""Tests binding prediction functions"""
	def setUp(self):
		""""""
		self.neoepitopes = {'CGCSQKCN': [('11', 71277056, '', 'AAA', 
										  'I', 0.1, 
										  'ENST00000398531.2_2')],
							'PVCCPCKI': [('11', 71277229, 'A', 'C', 'V', 
										  15.7, 'ENST00000398531.2_2')]}
		self.tools = {'mhcflurry1': ['mhcflurry-predict',
									 ['affinity', 'rank']]}
		self.alleles = ['HLA-A*02:01', 'HLA-B*07:02']
	def test_binding_scores(self):
		new_neoepitopes = gather_binding_scores(self.neoepitopes, 
												self.tools, 
												self.alleles)
		self.assertEqual(new_neoepitopes['CGCSQKCN'], [('11', 71277056, 
													'', 'AAA', 
													'I', 0.1, 
													'ENST00000398531.2_2',
													'31216.539180760206', 
													'98.32912499999998', 
													'34006.993380658656', 
													'80.02275'
												)])
		self.assertEqual(new_neoepitopes['PVCCPCKI'], [('11', 71277229, 
													'A', 'C', 'V', 
													15.7, 
													'ENST00000398531.2_2',
													'15595.009921030196', 
													'26.75187499999999', 
													'29813.963222395083', 
													'46.76187499999999'
												)])
		self.assertEqual(11, len(new_neoepitopes['PVCCPCKI'][0]))
		self.assertEqual(11, len(new_neoepitopes['CGCSQKCN'][0]))
class TestOutput(unittest.TestCase):
	"""Tests function to write output"""
	def setUp(self):
		"""Sets up paths and dictionaries"""
		self.base_dir = neoepiscope_dir
		self.out_file = os.path.join(self.base_dir, 'tests', 'neoepiscope.out'
						)
		self.correct_out = os.path.join(self.base_dir, 'tests', 'expected.neoepiscope.out'
						)
		self.tools = {'netMHCpan4': ['netMHCpan', ['rank', 'affinity']],
					  'netMHCIIpan3': ['netMHCIIpan', ['rank']]}
		self.HLA_alleles = ['HLA*A01:01', 'HLA*A02:01']
		self.neoepitopes = {'CGCSQKCN': [('11', 71277056, '',
										  'AAA', 'I', 0.1, 
										  'ENST00000398531.2_2',
										  5, 10000.0, 1, 5, 150, 4),
										 ('4', 300000, 'A',
										  'T', 'V', 10.2, 
										  'ENST00000398554.1_1',
										  0.5, 100, 0.5, 3, 150, 4)],
							'PVCCPCKI': [('11', 71277229, 'A', 'C', 'V', 
										  15.7, 'ENST00000398531.2_2',
										  10, 50.57, 1.2, 3, 10.1, 7),
										 ('11', 71277229, 'A', 'C', 'V', 
										  20.3, 'ENST00000398200.3_1',
										  10, 50.57, 1.2, 3, 10.1, 7)]}
	def testwrite(self):
		"""Tests that output file is written correctly"""
		write_results(self.out_file, self.HLA_alleles, self.neoepitopes,
					  self.tools)
		self.assertTrue(filecmp.cmp(self.out_file, self.correct_out))
	def tearDown(self):
		"""Removes test file"""
		os.remove(self.out_file)
unittest.main()
