#!/usr/bin/env python
"""
nexepi

Identifies neoepitopes from DNA-seq, VCF, GTF, and Bowtie index.
"""
import bisect
import argparse
import bowtie_index
import sys
import math
import string
import copy
import pickle
import copy
import os
import random
import re
import collections
import tempfile
import subprocess
import string
from transcript import Transcript, gtf_to_cds, cds_to_tree, get_transcripts_from_tree
from operator import itemgetter
from sortedcontainers import SortedDict
from intervaltree import Interval, IntervalTree
#import Hapcut2interpreter as hap

_revcomp_translation_table = string.maketrans('ATCG', 'TAGC')

# X below denotes a stop codon
_codon_table = {
        'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
        'TAT':'Y', 'TAC':'Y', 'TAA':'X', 'TAG':'X',
        'TGT':'C', 'TGC':'C', 'TGA':'X', 'TGG':'W',
        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
        'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
        'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }
_complement_table = string.maketrans('ATCG', 'TAGC')

_help_intro = '''neoepiscope searches for neoepitopes in seq data.'''
def help_formatter(prog):
    """ So formatter_class's max_help_position can be changed. """
    return argparse.HelpFormatter(prog, max_help_position=40)

def adjust_tumor_column(in_vcf, out_vcf):
    """ Swaps the sample columns in a somatic vcf

        HAPCUT2 only takes data from the first VCF sample column, so if the 
            tumor sample data is in the second VCF sample column, it must be
            swapped prior to optional germline merging or running HAPCUT2

        in_vcf: input vcf that needs the tumor sample data flipped
        out_vcf: output vcf to have the correct columns

        No return value.
    """
    header_lines = []
    other_lines = []
    # Process input vcf
    with open(in_vcf, 'r') as f:
        for line in f:
            # Preserve header lines with out change
            if line[0:2] == '##':
                header_lines.append(line.strip('\n'))
            # Adjust column header and variant lines
            else:
                tokens = line.strip('\n').split('\t')
                if line[0] == '#':
                    print ''.join(['WARNING: Reading ', tokens[9], 
                                    'as normal tissue and ', tokens[10],
                                    'as tumor tissue'])
                new_line = '\t'.join([tokens[0], tokens[1], tokens[2], 
                                        tokens[3], tokens[4], tokens[5], 
                                        tokens[6], tokens[7], tokens[8], 
                                        tokens[10], tokens[9]])
                other_lines.append(new_line)
    # Write new vcf
    with open(out_vcf, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        for line in other_lines:
            f.write(line + '\n')

def combinevcf(vcf1, vcf2, outfile='Combined.vcf'):
    """ Combines VCFs

        No return value.
    """
    vcffile = open(vcf2, 'r')
    temp = open(vcf2 + '.tumortemp', 'w+');
    header = open(vcf2 + '.header', 'w+');
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
        else:
            header.write(lines)
    vcffile.close()
    temp.close()
    header.close()
    vcffile = open(vcf1, 'r')
    temp = open(vcf2 + '.germlinetemp', 'w+');
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
    vcffile.close()
    temp.close()    
    markgermline = ''.join(['''awk '{print $0"*"}' ''', vcf2, 
                            ".germlinetemp > ", vcf2, '.germline'])
    marktumor    = ''.join(['''awk '{print $0}' ''', vcf2, 
                            '.tumortemp > ', vcf2, '.tumor'])
    subprocess.call(markgermline, shell=True)
    subprocess.call(marktumor, shell=True)
    command = ''.join(['cat ', vcf2, '.germline ', vcf2, '.tumor > ', 
                        vcf2, '.combine1'])
    subprocess.call(command, shell=True)
    command2 = ''.join(['sort -k1,1 -k2,2n ', vcf2, '.combine1 > ', 
                        vcf2, '.sorted'])
    subprocess.call(command2, shell=True)
    command3 = ''.join(['cat ', vcf2, '.header ', vcf2, '.sorted > ', 
                        vcf2, '.combine2'])
    subprocess.call(command3, shell=True)
    cut = ''.join(['cut -f1,2,3,4,5,6,7,8,9,10 ', vcf2, 
                    '.combine2 > ', outfile])
    subprocess.call(cut, shell=True)
    for file in ['.tumortemp', '.germlinetemp', '.combine1', '.combine2', 
                    '.sorted', '.tumor', '.germline', '.header']:
        cleanup = ''.join(['rm ', vcf2, file])
        subprocess.call(cleanup, shell=True)

def which(path):
    """ Searches for whether executable is present and returns version

        path: path to executable

        Return value: None if executable not found, else string with software
            name and version number
    """
    try:
        subprocess.check_call([path])
    except OSError as e:
        return None
    else:
        return path

def get_VAF_pos(VCF):
    """ Obtains position in VCF format/genotype fields of VAF

        VCF: path to input VCF

        Return value: None if VCF does not contain VAF, 
                        otherwise position of VAF
    """
    VAF_check = False
    with open(VCF) as f:
        for line in f:
            # Check header lines to see if FREQ exits in FORMAT fields
            if line[0] == '#':
                if 'FREQ' in line:
                    VAF_check = True
            else:
                # Check first entry to get position of FREQ if it exists
                if VAF_check:
                    tokens = line.strip('\n').split('\t')
                    format_field = tokens[8].split(':')
                    for i in range(0,len(format_field)):
                        if format_field[i] == 'FREQ':
                            VAF_pos = i
                            break
                # Return None if VCF does not contain VAF data
                else:
                    VAF_pos = None
                    break
    return VAF_pos

def process_haplotypes(hapcut_output, interval_dict):
    """ Stores all haplotypes relevant to different transcripts as a dictionary

        hapcut_output: output from HAPCUT2, adjusted to include unphased 
                        mutations as their own haplotypes (performed in 
                        software's prep mode)
        interval_dict: dictionary linking genomic intervals to transcripts

        Return value: dictinoary linking haplotypes to transcripts
    """
    affected_transcripts = collections.defaultdict(list)
    with open(hapcut_output, 'r') as f:
        block_transcripts = collections.defaultdict(list)
        for line in f:
            if line.startswith('BLOCK'):
                # Skip block header lines
                continue
            elif line[0] == '*':
                # Process all transcripts for the block
                for transcript_ID in block_transcripts:
                    block_transcripts[transcript_ID].sort(key=itemgetter(1))
                    for mut in block_transcripts[transcript_ID]:
                        affected_transcripts[transcript_ID].append(mut)
                # Reset transcript dictionary
                block_transcripts = collections.defaultdict(list)
            else:
                # Add mutation to transcript dictionary for the block
                tokens = line.strip("\n").split()
                if len(tokens[5]) == len(tokens[6]):
                    mutation_type = 'V'
                    pos = int(tokens[4])
                    ref = tokens[5]
                    alt = tokens[6]
                    mut_size = len(tokens[5])
                    end = pos + mut_size
                elif len(tokens[5]) > len(tokens[6]):
                    mutation_type = 'D'
                    deletion_size = len(tokens[5]) - len(tokens[6])
                    pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                    ref = tokens[5]
                    alt = deletion_size
                    end = pos + deletion_size
                elif len(tokens[5]) < len(tokens[6]):
                    mutation_type = 'I'
                    insertion_size = len(tokens[6]) - len(tokens[5])
                    pos = int(tokens[4])
                    ref = tokens[5]
                    alt = tokens[6][len(ref):]
                    end = pos + 1
                overlapping_transcripts = get_transcripts_from_tree(tokens[3], 
                                                                pos, 
                                                                end,
                                                                interval_dict)
                # For each overlapping transcript, add mutation entry
                # Contains chromosome, position, reference, alternate, allele
                #   A, allele B, genotype line from VCF
                for transcript in overlapping_transcripts:
                    block_transcripts[transcript].append([tokens[3], pos, 
                                                          ref, alt, 
                                                          tokens[1], tokens[2], 
                                                          tokens[7], mutation_type])
    return affected_transcripts

def get_affinity_netMHCIIpan(peptides, allele, netmhciipan=program,
                                            remove_files=True):
                """ Obtains binding affinities from list of peptides

                    peptides: peptides of interest (list of strings)
                    allele: Allele to use for binding affinity (string)
                    remove_files: option to remove intermediate files

                    Return value: affinities (a list of binding affinities 
                                    as strings)
                """
                files_to_remove = []
                try:
                    # Check that allele is valid for method
                    with open(''.join([os.path.dirname(__file__),
                             '/availableAlleles.pickle']), 'rb'
                            ) as allele_stream:
                        avail_alleles = pickle.load(allele_stream)
                    # Homogenize format
                    allele = allele.replace('HLA-', '')
                    if allele not in avail_alleles['netMHCIIpan']:
                        print ''.join(['WARNING: ', allele, 
                                    ' is not a valid allele for netMHCIIpan'])
                        return ['NA' for i in range(0,len(peptides))]
                    # Establish return list and sample id
                    sample_id = '.'.join([peptides[0],
                                            str(len(peptides)), allele,
                                            'netmhciipan'])
                    affinities = []
                    # Write one peptide per line to a temporary file for 
                    #   input if peptide length is at least 9
                    # Count instances of smaller peptides
                    na_count = 0
                    peptide_file = tempfile.mkstemp(
                                    suffix='.peptides', prefix='id.', text=True)
                    files_to_remove.append(peptide_file)
                    with open(peptide_file[1], 'w') as f:
                        for sequence in peptides:
                            if len(sequence) >= 9:
                                print >>f, sequence
                            else:
                                na_count += 1
                    if na_count > 0:
                        print ' ' .join(['WARNING:', str(na_count),
                                        'peptides not compatible with'
                                        'netMHCIIpan will not receive score'])
                    # Establish temporary file to hold output
                    mhc_out = tempfile.mkstemp(suffix='.netMHCIIpan.out', 
                                                prefix='id.', text=True)
                    files_to_remove.append(mhc_out)
                    # Run netMHCIIpan
                    subprocess.check_call(
                                    [netmhciipan, '-a', allele, '-inptype', '1', 
                                     '-xls', '-xlsfile', mhc_out, peptide_file]
                                )
                    # Retrieve scores for valid peptides
                    score_dict = {}
                    with open(mhc_out[1], 'r') as f:
                        # Skip headers
                        f.readline()
                        f.readline()
                        for line in f:
                            tokens = line.split('\t')
                            # token 1 is peptide; token 4 is score
                            score_dict[tokens[1]] = tokens[4]
                    # Produce list of scores for valid peptides
                    # Invalid peptides receive "NA" score
                    for sequence in peptides:
                        if sequence in score_dict:
                            nM = score_dict[sequence]
                        else:
                            nM = 'NA'
                            affinities.append(nM)
                    return affinities
                finally:
                    if remove_files:
                        for file_to_remove in files_to_remove:
                            os.remove(file_to_remove)

def get_affinity_netMHCpan(peptides, allele, netmhcpan=program,
                                remove_files=True):
                """ Obtains binding affinities from list of peptides

                    peptides: peptides of interest (list of strings)
                    allele: allele to use for binding affinity 
                                (string, format HLA-A02:01)
                    remove_files: option to remove intermediate files

                    Return value: affinities (a list of binding affinities 
                                    as strings)
                """
                files_to_remove = []
                try:
                    # Check that allele is valid for method
                    with open(''.join([os.path.dirname(__file__),
                             '/availableAlleles.pickle']), 'rb'
                            ) as allele_stream:
                        avail_alleles = pickle.load(allele_stream)
                    allele = allele.replace('*', '')
                    if allele not in avail_alleles['netMHCpan']:
                        print ''.join(['WARNING: ', allele, 
                                    ' is not a valid allele for netMHCpan'])
                        return ['NA' for i in range(0,len(peptides))]
                    # Establish return list and sample id
                    sample_id = '.'.join([peptides[0], str(len(peptides)), 
                                            allele, 'netmhcpan'])
                    affinities = []
                    # Write one peptide per line to a temporary file for input
                    peptide_file = tempfile.mkstemp(suffix='.peptides', 
                                                    prefix=''.join([sample_id, 
                                                                    '.']), 
                                                    text=True)
                    files_to_remove.append(peptide_file)
                    with open(peptide_file[1], 'w') as f:
                        for sequence in peptides:
                            print >>f, sequence
                    # Establish temporary file to hold output
                    mhc_out = tempfile.mkstemp(suffix='.netMHCpan.out', 
                                                prefix=''.join([sample_id, 
                                                                '.']), 
                                                text=True)
                    files_to_remove.append(mhc_out)
                    # Run netMHCpan
                    subprocess.check_call(
                        [netmhcpan, '-a', allele, '-inptype', '1', '-p', '-xls', 
                            '-xlsfile', mhc_out, peptide_file])
                    with open(mhc_out[1], 'r') as f:
                        f.readline()
                        f.readline()
                        for line in f:
                            line = line.strip('\n').split('\t')
                            nM = line[5]
                            affinities.append(nM)
                    return affinities
                finally:
                    # Remove temporary files
                    if remove_files:
                        for file_to_remove in files_to_remove:
                            os.remove(file_to_remove)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro, 
                                        formatter_class=help_formatter)
    subparsers = parser.add_subparsers(help=(
                                    'subcommands; add "-h" or "--help" '
                                    'after a subcommand for its parameters'
                                ), dest='subparser_name')
    test_parser = subparsers.add_parser('test', help='performs unit tests')
    index_parser = subparsers.add_parser('index',
                                        help=('produces pickled dictionaries '
                                        'linking transcripts to intervals and '
                                        ' CDS lines in a GTF'))
    swap_parser = subparsers.add_parser('swap',
                                        help=('swaps tumor and normal columns '
                                        'in a somatic vcf if necessary for '
                                        'proper HapCUT2 results'))
    merge_parser = subparsers.add_parser('merge',
                                         help=('merges germline and somatic '
                                               'VCFS for combined mutation '
                                               'phasing with HAPCUT2'))
    prep_parser = subparsers.add_parser('prep',
                                        help=('combines HAPCUT2 output with '
                                              'unphased variants for call '
                                              'mode'))
    call_parser = subparsers.add_parser('call', help='calls neoepitopes')
    # Index parser options (produces pickled dictionaries for transcript data)
    index_parser.add_argument('-g', '--gtf', type=str, required=True,
            help='input path to GTF file'
        )  
    index_parser.add_argument('-d', '--dicts', type=str, required=True,
            help='output path to pickled CDS dictionary'
        )
    # Swap parser options (swaps columns in somatic VCF)
    swap_parser.add_argument('-i', '--input', type=str, required=True,
            help='input path to somatic VCF'
        )
    swap_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to column-swapped VCF'
        )
    # Merger parser options (merges somatic and germline VCFs)
    merge_parser.add_argument('-g', '--germline', type=str, required=True,
            help='input path to germline VCF'
        )
    merge_parser.add_argument('-s', '--somatic', type=str, required=True,
            help='input path to somatic VCF'
        )
    merge_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to combined VCF'
        )
    # Prep parser options (adds unphased mutations as their own haplotype)
    prep_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input VCF'
        )
    prep_parser.add_argument('-c', '--hapcut2-output', type=str, required=True,
            help='path to output file of HAPCUT2 run on input VCF'
        )
    prep_parser.add_argument('-o', '--output', type=str, required=True,
            help='path to output file to be input to call mode'
        )
    # Call parser optinos (calls neoepitopes)
    call_parser.add_argument('-x', '--bowtie-index', type=str, required=True,
            help='path to Bowtie index basename'
        )
    call_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input path to VCF'
        )
    call_parser.add_argument('-d', '--dicts', type=str, required=True,
            help='input path to pickled CDS dictionary'
        )
    call_parser.add_argument('-c', '--merged-hapcut2-output', type=str,
            required=True,
            help='path to output of prep subcommand'
        )
    call_parser.add_argument('-k', '--kmer-size', type=str, required=False,
            default='8,11', help='kmer size for epitope calculation'
        )
    call_parser.add_argument('-p', '--affinity-predictor', type=str, 
            required=False, default='netMHCpan', 
            help='path to executable for binding affinity prediction software; '
                'for multiple softwares, use comma-separated list of paths'
        )
    call_parser.add_argument('-s', '--predictor-scores', type=str, 
            required=False, default='rank', 
            help='comma separated list of scoring methods to use, see '
            'documentation online for more information'
        )
    call_parser.add_argument('-a', '--alleles', type=str, required=True,
            help='comma separated list of alleles; '
                 'see documentation online for more information'
        )
    call_parser.add_argument('-o', '--output_file', type=str, required=True,
            help='path to output file'
        )
    call_parser.add_argument('-u', '--upstream_atgs', type=str, required=False,
            default='novel', help='how to handle upstream start codons, see '
            'documentation online for more information'
        )
    args = parser.parse_args()
    
    if args.subparser_name == 'test':
        del sys.argv[1:] # Don't choke on extra command-line parameters
        import unittest
        import filecmp
        class TestGTFprocessing(unittest.TestCase):
            """Tests proper creation of dictionaries store GTF data"""
            def setUp(self):
                """Sets up gtf file and creates dictionaries for tests"""
                self.gtf = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.gtf'
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
                self.varscan = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.varscan.vcf'
                                )                
                self.germline = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.germline.vcf'
                                )   
                self.precombined = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.combined.vcf'
                                )  
                self.outvcf = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.testcombine.vcf'
                                )  
                combinevcf(self.germline, self.varscan, self.outvcf)
            def test_merge(self):
                """Fails if VCFs were merged improperly"""
                self.assertTrue(filecmp.cmp(self.outvcf, self.precombined))
            def tearDown(self):
                """Removes test file"""
                os.remove(self.outvcf)
        class TestVAFpos(unittest.TestCase):
            """Tests fetching of VAF position from VCF file"""
            def setUp(self):
                """ Sets up vcf files to use for tests """
                self.varscan = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.varscan.vcf'
                                )
                self.mutect = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.mutect.vcf'
                                )
            def test_position(self):
                """Fails if incorrect positions are returned"""
                self.assertEqual(get_VAF_pos(self.varscan), 5)
                self.assertEqual(get_VAF_pos(self.mutect), None)
        class TestHaplotypeProcessing(unittest.TestCase):
            """Tests proper processing of HAPCUT2 files"""
            def setUp(self):
                """Sets up input files and dictionaries to use for tests"""
                self.Ygtf = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.gtf'
                                )
                self.Ycds = gtf_to_cds(self.Ygtf, 'NA', pickle_it=False)
                self.Ytree = cds_to_tree(self.Ycds, 'NA', pickle_it=False)
                self.Yhapcut = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.hapcut.out'
                                )
                self.Chr11gtf = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Chr11.gtf'
                                )
                self.Chr11cds = gtf_to_cds(self.Chr11gtf, 'NA', pickle_it=False)
                self.Chr11tree = cds_to_tree(self.Chr11cds, 'NA',
                                                pickle_it=False)
                self.Chr11hapcut = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Chr11.hapcut.out'
                                ) 
            def test_hap_processing(self):
                """Fails if file is processed incorrectly"""
                Chr11_txs = process_haplotypes(self.Chr11hapcut, self.Chr11tree)
                self.assertEqual(sorted(Chr11_txs.keys()), 
                                        ['ENST00000398531.2_2'])
                self.assertEqual(Chr11_txs['ENST00000398531.2_2'],
                                [['11', 71276862, 'TGT', 2, '0', '0', 
                                   '0/0:.:53:52:0:0%:22,30,0,0:.:2', 'D'], 
                                 ['11', 71276900, 'C', 'G', '0', '0',
                                  '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'V'],
                                  ['11', 71277000, 'G', 'AA', '0', '0',
                                  '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'I']]
                    )
        unittest.main()
    elif args.subparser_name == 'index':
        cds_dict = gtf_to_cds(args.gtf, args.dicts)
        tree = cds_to_tree(cds_dict, args.dicts)
        # FM indexing of proteome??
    elif args.subparser_name == 'swap':
        adjust_tumor_column(args.input, args.output)
    elif args.subparser_name == 'merge':
        combinevcf(args.germline, args.somatic, outfile=args.output)
    elif args.subparser_name == 'prep':
        phased = collections.defaultdict(set)
        with open(args.output, 'w') as output_stream:
            print >>output_stream, '********'
            with open(args.hapcut2_output) as hapcut2_stream:
                for line in hapcut2_stream:
                    if line[0] != '*' and not line.startswith('BLOCK'):
                        tokens = line.strip().split('\t')
                        phased[(tokens[3], int(tokens[4]))].add(
                                                        (tokens[5], tokens[6])
                                                    )
                    print >>output_stream, line,
            with open(args.vcf) as vcf_stream:
                first_char = '#'
                while first_char == '#':
                    line = vcf_stream.readline().strip()
                    try:
                        first_char = line[0]
                    except IndexError:
                        first_char = '#'
                counter = 1
                while line:
                    tokens = line.split('\t')
                    pos = int(tokens[1])
                    alt_alleles = tokens[4].split(',')
                    for allele in alt_alleles:
                        if (tokens[3], allele) not in phased[
                                                    (tokens[0], pos)
                                                ]:
                            print >>output_stream, 'BLOCK: unphased'
                            print >>output_stream, ('{vcf_line}\tNA\tNA\t{chrom}\t'
                                                   '{pos}\t{ref}\t{alt}\t'
                                                   '{genotype}\tNA\tNA\tNA').format(
                                                        vcf_line=counter,
                                                        chrom=tokens[0],
                                                        pos=pos,
                                                        ref=tokens[3],
                                                        alt=tokens[4],
                                                        genotype=tokens[9]
                                                    )
                            print >>output_stream, '********' 
                    line = vcf_stream.readline().strip()
                    counter += 1
    elif args.subparser_name == 'call':
        # Load pickled dictionaries
        interval_dict = pickle.load(open(args.dicts + ''.join([dictdir, 
                                '/', 'intervals_to_transcript.pickle']), 'rb'))
        cds_dict = pickle.load(open(args.dicts + ''.join([dictdir, 
                                '/', 'transcript_to_CDS.pickle']), 'rb'))
        # Check affinity predictor
        for tool in ','.split(args.affinity_predictor):
            program = which(tool)
            if program is None:
                raise ValueError(' '.join([program, 'is not a valid software']))
            elif 'netMHCIIpan' in program:
                pass
                # Deal with version/score checking
            elif 'netMHCpan' in program:
                pass
                # Deal with version/score checking
            else:
                raise ValueError(' '.join([program, 'is not a valid software']))
        # Obtain VAF frequency VCF position
        VAF_pos = get_VAF_pos(args.vcf)
        # Obtain peptide sizes for kmerizing peptides
        if ',' in args.kmer_size:
            size_list = args.kmer_size.split(',')
            size_list.sort()
            for i in range(0, len(size_list)):
                size_list[i] = int(size_list[i])
        hla_alleles = args.alleles.split(',')
        # For retrieving genome sequence
        reference_index = bowtie_index.BowtieIndexReference(args.bowtie_index)
        # Find transcripts that haplotypes overlap 
        relevant_transcripts = process_haplotypes(args.merged_hapcut2_output, 
                                                    interval_dict)
        # Iterate over relevant transcripts to create transcript objects and
        #   enumerate neoepitopes
        neoepitopes = collections.defaultdict(list)
        # Establish handling of ATGs
        if arg.upstream_atgs == 'novel':
            only_novel_upstream = True
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == 'all':
            only_novel_upstream = False
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == 'none':
            only_novel_upstream = False
            only_downstream = True
            only_reference = False
        elif args.upstream_atgs == 'reference':
            only_novel_upstream = False
            only_downstream = False
            only_reference = True
        else:
            raise RuntimeError('--upstream_atgs must be one of {"novel", "all", '
                                '"none", "reference"}')
        for affected_transcript in relevant_transcripts:
            # Create transcript object
            transcriptA = Transcript(reference_index, 
                            [[str(chrom), 'blah', seq_type, str(start), 
                              str(end), '.', strand] for (chrom, seq_type, 
                                                          start, end, strand) 
                              in cds_dict[transcript_ID]], transcript_ID
                            )
            transcriptB = Transcript(reference_index, 
                            [[str(chrom), 'blah', seq_type, str(start), 
                              str(end), '.', strand] for (chrom, seq_type, 
                                                          start, end, strand) 
                              in cds_dict[transcript_ID]], transcript_ID
                            )
            # Iterate over haplotypes associated with this transcript
            haplotypes = relevant_transcripts[affected_transcript]
            for ht in haplotypes:
                # Make edits for each mutation
                for mutation in ht:
                    # Determine if mutation is somatic or germline
                    if mutation[6][-1] == '*':
                        mutation_class = 'G'
                    else:
                        mutation_class = 'S'
                    # Determine VAF if available
                    if VAF_pos is not None:
                        VAF = float(mutation[6][VAF_pos].strip('*').strip('%'))
                    else:
                        VAF = None
                    # Determine which copies variant exists on & make edits
                    if mutation[4] == '1':
                        transcriptA.edit(mutation[3], mutation[1], 
                                    mutation_type=mutation[7], 
                                    mutation_class=mutation_class,
                                    VAF=VAF)
                    if mutation[5] == '1':
                        transcriptB.edit(mutation[3], mutation[1], 
                                    mutation_type=mutation[7], 
                                    mutation_class=mutation_class,
                                    VAF=VAF)
                # Extract neoepitopes
                A_peptides = transcriptA.neopeptides(min_size=size_list[0], 
                                                     max_size=size_list[-1],
                                                     include_somatic=1,
                                                     include_germline=2, 
                                                     only_novel_upstream=only_novel_upstream,
                                                     only_downstream=only_downstream, 
                                                     only_reference=only_reference)
                B_peptides = transcriptB.neopeptides(min_size=size_list[0], 
                                                     max_size=size_list[-1],
                                                     include_somatic=1,
                                                     include_germline=2, 
                                                     only_novel_upstream=only_novel_upstream,
                                                     only_downstream=only_downstream, 
                                                     only_reference=only_reference)
                # Store neoepitopes and their metadata
                for pep in A_peptides:
                    for meta_data in A_peptides[pep]:
                        if meta_data not in neoepitopes[pep]:
                            neoepitopes[pep].append(meta_data + (transcriptA.transcript_ID,))
                for pep in B_peptides:
                    for meta_data in B_peptides[pep]:
                        if meta_data not in neoepitopes[pep]:
                            neoepitopes[pep].append(meta_data + (transcriptB.transcript_ID,))
                transcriptA.reset(reference=True)
                transcriptB.reset(reference=True)
        for allele in hla_alleles:
            ## What value do we want from netMHCpan? affinity or score?
            binding_scores = get_affinity(sorted(neoepitopes.keys()), allele)
            for i in range(0, sorted(neoepitopes.keys())):
                meta_data = neoepitopes[sorted(neoepitopes.keys())[i]]
                for mutation in meta_data:
                    mutation = mutation + (binding_scores[i],)
        ## What format do we want for the output file?
        with open(args.output_file, 'w') as o:
            headers = ['Neoepitope', 'Chromsome', 'Pos', 'Ref', 'Alt', 'VAF', 
                        'Transcript_ID']
            for allele in hla_alleles:
                headers.append('_'.join([allele, 'score']))
            o.write('\t'.join(headers) + '\n')
            for epitope in neoepitopes:
                for mutation in neoepitopes[epitope][2]:
                    out_line = [epitope, mutation[1], str(mutation[2]), mutation[3], 
                                neoepitopes[epitope][0], str(mutation[4]),
                                mutation[5]]
                    for i in range(5:len(mutation)):
                        out_line.append(str(mutation[i]))
                o.write('\t'.join(out_line) + '\n')
    else:
        sys.exit(''.join([args.subparser_name, 
                            ' is not a valid software mode']))
