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
        "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
        "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
        "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
        "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }
_complement_table = string.maketrans("ATCG", "TAGC")

_help_intro = """neoscan searches for neoepitopes in seq data."""
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

def combinevcf(vcf1, vcf2, outfile="Combined.vcf"):
    """ Combines VCFs

        No return value.
    """
    vcffile = open(vcf2, "r")
    temp = open(vcf2 + ".tumortemp", "w+");
    header = open(vcf2 + ".header", "w+");
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
        else:
            header.write(lines)
    vcffile.close()
    temp.close()
    header.close()
    vcffile = open(vcf1, "r")
    temp = open(vcf2 + ".germlinetemp", "w+");
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
    vcffile.close()
    temp.close()    
    markgermline = "".join(['''awk '{print $0"*"}' ''', vcf2, 
                            ".germlinetemp > ", vcf2, ".germline"])
    marktumor    = "".join(['''awk '{print $0}' ''', vcf2, 
                            ".tumortemp > ", vcf2, ".tumor"])
    subprocess.call(markgermline, shell=True)
    subprocess.call(marktumor, shell=True)
    command = "".join(["cat ", vcf2, ".germline ", vcf2, ".tumor > ", 
                        vcf2, ".combine1"])
    subprocess.call(command, shell=True)
    command2 = "".join(["sort -k1,1 -k2,2n ", vcf2, ".combine1 > ", 
                        vcf2, ".sorted"])
    subprocess.call(command2, shell=True)
    command3 = "".join(["cat ", vcf2, ".header ", vcf2, ".sorted > ", 
                        vcf2, ".combine2"])
    subprocess.call(command3, shell=True)
    cut = "".join(["cut -f1,2,3,4,5,6,7,8,9,10 ", vcf2, 
                    ".combine2 > ", outfile])
    subprocess.call(cut, shell=True)
    for file in [".tumortemp", ".germlinetemp", ".combine1", ".combine2", 
                    ".sorted", ".tumor", ".germline", ".header"]:
        cleanup = "".join(["rm ", vcf2, file])
        subprocess.call(cleanup,shell=True)

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
            if line[0] == "#":
                if "FREQ" in line:
                    VAF_check = True
            else:
                # Check first entry to get position of FREQ if it exists
                if VAF_check:
                    tokens = line.strip("\n").split("\t")
                    format_field = tokens[8].split(":")
                    for i in range(0,len(format_field)):
                        if format_field[i] == "FREQ":
                            VAF_pos = i
                            break
                # Return None if VCF does not contain VAF data
                else:
                    VAF_pos = None
                    break
    return VAF_pos

def seq_to_peptide(seq, reverse_strand=False):
    """ Translates nucleotide sequence into peptide sequence.

        All codons after stop codon (X) are omitted.

        seq: nucleotide sequence
        reverse_strand: True iff strand is -

        Return value: peptide string
    """
    seq_size = len(seq)
    if reverse_strand:
        seq = seq[::-1].translate(_complement_table)
    peptide = []
    for i in xrange(0, seq_size - seq_size % 3, 3):
        codon = _codon_table[seq[i:i+3]]
        peptide.append(codon)
        if codon == 'X':
            break
    # defunct behavior: return X's for all codons after stop codon
    #for j in xrange(i + 3, seq_size - seq_size % 3, 3):
    #    peptide.append('X')
    return ''.join(peptide)

def median(numbers):
    numbers = sorted(numbers)
    center = len(numbers) / 2
    if len(numbers) == 1:
        return numbers[0]
    elif len(numbers) % 2 == 0:
        return sum(numbers[center - 1:center + 1]) / 2.0
    else:
        return numbers[center]

def kmerize_peptide(peptide, min_size=8, max_size=11):
    """ Obtains subsequences of a peptide.

        normal_peptide: normal peptide seq
        min_size: minimum subsequence size
        max_size: maximum subsequence size

        Return value: list of all possible subsequences of size between
            min_size and max_size
    """
    peptide_size = len(peptide)
    return [item for sublist in
                [[peptide[i:i+size] for i in xrange(peptide_size - size + 1)]
                    for size in xrange(min_size, max_size + 1)]
            for item in sublist if 'X' not in item]

def neoepitopes(mutated_seq, reverse_strand=False, min_size=8, max_size=11):
    """ Finds neoepitopes from normal and mutated seqs.

        mutated_seq: mutated nucelotide sequence
        reverse_strand: True iff strand is -
        min_size: minimum peptide kmer size to write
        max_size: maximum petide kmer size to write

        Return value: List of mutated_kmers
    """
    return kmerize_peptide(seq_to_peptide(mutated_seq, 
                                            reverse_strand=reverse_strand
                                            ),
                            min_size=min_size,
                            max_size=max_size
                            )

def process_haplotypes(hapcut_output, interval_dict):
    """ Stores all haplotypes relevant to different transcripts as a dictionary

        hapcut_output: output from HAPCUT2, adjusted to include unphased 
                        mutations as their own haplotypes (performed in 
                        software's prep mode)
        interval_dict: dictionary linking genomic intervals to transcripts

        Return value: dictinoary linking haplotypes to transcripts
    """
    affected_transcripts = collections.defaultdict(list)
    with open(hapcut_output, "r") as f:
        block_transcripts = collections.defaultdict(list)
        for line in f:
            if line.startswith('BLOCK'):
                # Skip block header lines
                continue
            elif line[0] == "*":
                # Process all transcripts for the block
                for transcript_ID in block_transcripts:
                    block_transcripts[transcript_ID].sort(key=itemgetter(1))
                    affected_transcripts[transcript_ID].append(
                                                block_transcripts[transcript_ID]
                                                )
                # Reset transcript dictionary
                block_transcripts = collections.defaultdict(list)
            else:
                # Add mutation to transcript dictionary for the block
                tokens = line.strip("\n").split()
                mut_size = min(len(tokens[5]), len(tokens[6]))
                end = int(tokens[4]) + mut_size
                overlapping_transcripts = get_transcripts_from_tree(tokens[3], 
                                                                int(tokens[4]), 
                                                                end,
                                                                interval_dict)
                # For each overlapping transcript, add mutation entry
                # Contains chromosome, position, reference, alternate, allele
                #   A, allele B, genotype line from VCF
                for transcript in overlapping_transcripts:
                    block_transcripts[transcript].append([tokens[3], tokens[4], 
                                                          tokens[5], tokens[6], 
                                                          tokens[1], tokens[2], 
                                                          tokens[7]])
    return affected_transcripts


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
    index_parser.add_argument('-g', '--gtf', type=str, required=True,
            help='input path to GTF file'
        )  
    index_parser.add_argument('-d', '--dicts', type=str, required=True,
            help='output path to pickled CDS dictionary'
        )
    swap_parser.add_argument('-i', '--input', type=str, required=True,
            help='input path to somatic VCF'
        )
    swap_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to column-swapped VCF'
        )
    merge_parser.add_argument('-g', '--germline', type=str, required=True,
            help='input path to germline VCF'
        )
    merge_parser.add_argument('-s', '--somatic', type=str, required=True,
            help='input path to somatic VCF'
        )
    merge_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to combined VCF'
        )
    prep_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input VCF'
        )
    prep_parser.add_argument('-c', '--hapcut2-output', type=str, required=True,
            help='path to output file of HAPCUT2 run on input VCF'
        )
    prep_parser.add_argument('-o', '--output', type=str, required=True,
            help='path to output file to be input to call mode'
        )
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
    call_parser.add_argument('-m', '--method', type=str, required=False,
            default='-', 
            help='method for calculating epitope binding affinities'
        )
    call_parser.add_argument('-p', '--affinity-predictor', type=str, 
            required=False, default='netMHCpan', 
            help='path to executable for binding affinity prediction software'
        )
    call_parser.add_argument('-a', '--allele', type=str, required=True,
            help='allele; see documentation online for more information'
        )
    call_parser.add_argument('-f', '--VAF-freq-calc', type=str, required=False,
            default='median',
            help='method for calculating VAF: choice of mean, median, min, max'
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
                self.assertEqual(len(self.Ycds.keys()), 901)
            def test_CDS_tree(self):
                """Fails if dictionary was built incorrectly"""
                self.assertEqual(len(self.Ytree.keys()), 1)
                self.assertEqual(len(self.Ytree['Y']), 4808)
            def test_transcript_extraction(self):
                """Fails if incorrect transcripts are pulled"""
                self.assertEqual(len(get_transcripts_from_tree(
                                                                'Y', 150860, 
                                                                150861,
                                                                self.Ytree)), 
                                                                11
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
                                                    'ENST00000399012.6_3_PAR_Y',
                                                    'ENST00000415337.6_3_PAR_Y',
                                                    'ENST00000429181.6_2_PAR_Y',
                                                    'ENST00000430923.7_3_PAR_Y',
                                                    'ENST00000443019.6_2_PAR_Y',
                                                    'ENST00000445062.6_2_PAR_Y',
                                                    'ENST00000447472.6_3_PAR_Y',
                                                    'ENST00000448477.6_2_PAR_Y',
                                                    'ENST00000484611.7_1_PAR_Y']
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
                Ychrom_txs = process_haplotypes(self.Yhapcut, self.Ytree)
                Chr11_txs = process_haplotypes(self.Chr11hapcut, self.Chr11tree)
                self.assertEqual(sorted(Ychrom_txs.keys()), 
                                        ['ENST00000431853.1_1'])
                self.assertEqual(Ychrom_txs['ENST00000431853.1_1'], 
                                 [[['Y', '59001513', 'G', 'A', '0', '0', 
                                    '0/0:.:256:242:13:5.1%:184,58,9,4:.:2'], 
                                 ['Y', '59001559', 'G', 'A', '0', '0',
                                  '0/0:.:276:261:14:5.09%:117,144,6,8:.:2']]]
                                )
                self.assertEqual(sorted(Chr11_txs.keys()), 
                                        ['ENST00000398531.2_2', 
                                        'ENST00000528813.1_1'])
                self.assertEqual(Chr11_txs['ENST00000398531.2_2'],
                                [[['11', '71276861', 'T', 'C', '0', '0', 
                                   '0/0:.:53:52:0:0%:22,30,0,0:.:2'], 
                                 ['11', '71276900', 'C', 'G', '0', '0',
                                  '0/0:.:35:34:0:0%:19,15,0,0:.:2']]]
                    )
                self.assertEqual(Chr11_txs['ENST00000528813.1_1'],
                                [[['11', '56099004', 'C', 'T', '0', '0', 
                                   '0/0:.:30:30:0:0%:14,16,0,0:.:2'], 
                                 ['11', '56099010', 'A', 'G', '0', '0',
                                  '0/0:.:29:29:0:0%:15,14,0,0:.:2']]]
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
        interval_dict = pickle.load(open(args.dicts + "".join([dictdir, 
                                "/", "intervals_to_transcript.pickle"]), "rb"))
        cds_dict = pickle.load(open(args.dicts + "".join([dictdir, 
                                "/", "transcript_to_CDS.pickle"]), "rb"))
        # Check affinity predictor
        program = which(args.affinity_predictor)
        if program is None:
            raise ValueError(" ".join([program, "is not a valid software"]))
        elif "netMHCIIpan" in program:
            def get_affinity(peptides, allele, netmhciipan=program,
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
                    with open("".join([os.path.dirname(__file__),
                             "/availableAlleles.pickle"]), "rb"
                            ) as allele_stream:
                        avail_alleles = pickle.load(allele_stream)
                    # Homogenize format
                    allele = allele.replace("HLA-", "")
                    if allele not in avail_alleles["netMHCIIpan"]:
                        sys.exit(" ".join([allele,
                                 " is not a valid allele for netMHCIIpan"])
                                )
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
                                    suffix=".peptides", prefix="id.", text=True)
                    files_to_remove.append(peptide_file)
                    with open(peptide_file[1], "w") as f:
                        for sequence in peptides:
                            if len(sequence) >= 9:
                                print >>f, sequence
                            else:
                                na_count += 1
                    if na_count > 0:
                        print ' ' .join(['Warning:', str(na_count),
                                        'peptides not compatible with'
                                        'netMHCIIpan will not receive score'])
                    # Establish temporary file to hold output
                    mhc_out = tempfile.mkstemp(suffix=".netMHCIIpan.out", 
                                                prefix="id.", text=True)
                    files_to_remove.append(mhc_out)
                    # Run netMHCIIpan
                    subprocess.check_call(
                                    [netmhciipan, "-a", allele, "-inptype", "1", 
                                     "-xls", "-xlsfile", mhc_out, peptide_file]
                                )
                    # Retrieve scores for valid peptides
                    score_dict = {}
                    with open(mhc_out[1], "r") as f:
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
                            nM = "NA"
                            affinities.append(nM)
                    return affinities
                finally:
                    if remove_files:
                        for file_to_remove in files_to_remove:
                            os.remove(file_to_remove)
        elif "netMHCpan" in program:
            # define different affinity prediction function here
            def get_affinity(peptides, allele, netmhcpan=program,
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
                    with open("".join([os.path.dirname(__file__),
                             "/availableAlleles.pickle"]), "rb"
                            ) as allele_stream:
                        avail_alleles = pickle.load(allele_stream)
                    allele = allele.replace("*", "")
                    if allele not in avail_alleles["netMHCpan"]:
                        sys.exit(" ".join([allele,
                                 " is not a valid allele for netMHC"])
                                )
                    # Establish return list and sample id
                    sample_id = '.'.join([peptides[0], str(len(peptides)), 
                                            allele, 'netmhcpan'])
                    affinities = []

                    # Write one peptide per line to a temporary file for input
                    peptide_file = tempfile.mkstemp(suffix=".peptides", 
                                                    prefix="".join([sample_id, 
                                                                    "."]), 
                                                    text=True)
                    files_to_remove.append(peptide_file)
                    with open(peptide_file[1], "w") as f:
                        for sequence in peptides:
                            print >>f, sequence

                    # Establish temporary file to hold output
                    mhc_out = tempfile.mkstemp(suffix=".netMHCpan.out", 
                                                prefix="".join([sample_id, 
                                                                "."]), 
                                                text=True)
                    files_to_remove.append(mhc_out)
                    # Run netMHCpan
                    subprocess.check_call(
                        [netmhcpan, "-a", allele, "-inptype", "1", "-p", "-xls", 
                            "-xlsfile", mhc_out, peptide_file])
                    with open(mhc_out[1], "r") as f:
                        f.readline()
                        f.readline()
                        for line in f:
                            line = line.strip("\n").split("\t")
                            nM = line[5]
                            affinities.append(nM)
                    return affinities
                finally:
                    # Remove temporary files
                    if remove_files:
                        for file_to_remove in files_to_remove:
                            os.remove(file_to_remove)
        else:
            raise ValueError(" ".join([program, "is not a valid software"]))
        # Obtain VAF frequency VCF position
        VAF_pos = get_VAF_pos(args.vcf)
        # Obtain peptide sizes for kmerizing peptides
        if ',' in args.kmer_size:
            size_list = args.kmer_size.split(',')
            size_list.sort()
            for i in range(0, len(size_list)):
                size_list[i] = int(size_list[i])
        # For retrieving genome sequence
        reference_index = bowtie_index.BowtieIndexReference(args.bowtie_index)
        # Find transcripts that haplotypes overlap 
        relevant_transcripts = process_haplotypes(args.merged_hapcut2_output, 
                                                    interval_dict)
        # Iterate over relevant transcripts to create transcript objects and
        #   enumerate neoepitopes
        neoepitopes = collections.defaultdict(list)
        for affected_transcript in relevant_transcripts:
            # Create transcript object
            transcriptA = Transcript(reference_index, 
                            [[str(chrom), 'blah', seq_type, str(start), 
                              str(end), '.', strand] for (chrom, seq_type, 
                                                          start, end, strand) 
                              in cds_dict[transcript_ID]]
                            )
            transcriptB = Transcript(reference_index, 
                            [[str(chrom), 'blah', seq_type, str(start), 
                              str(end), '.', strand] for (chrom, seq_type, 
                                                          start, end, strand) 
                              in cds_dict[transcript_ID]]
                            )
            # Iterate over haplotypes associated with this transcript
            haplotypes = relevant_transcripts[affected_transcript]
            for ht in haplotypes:
                # Make edits for each mutation
                for mutation in ht:
                    # Determine type of mutation
                    if len(mutation[2] == len(mutation[3])):
                        mutation_type = 'V'
                    elif len(mutation[2]) < len(mutation[3]):
                        mutation_type = 'I'
                    elif len(mutation[2]) > len(mutation[3]):
                        mutation_type = 'D'
                    else:
                        mutation_type = '?'
                    # Determine if mutation is somatic or germline
                    if mutation[6][-1] == "*":
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
                        transciptA.edit(mutation[3], int(mutation[1]), 
                                    mutation_type=mutation_type, 
                                    mutation_class=mutation_class,
                                    VAF=VAF)
                    if mutation[5] == '1':
                        transciptB.edit(mutation[3], int(mutation[1]), 
                                    mutation_type=mutation_type, 
                                    mutation_class=mutation_class,
                                    VAF=VAF)
                # Extract neoepitopes
                A_peptides = transcriptA.neopeptides(min_size=size_list[0], 
                                                     max_size=size_list[-1],
                                                     include_somatic=1,
                                                     include_germline=2)
                B_peptides = transcriptB.neopeptides(min_size=size_list[0], 
                                                     max_size=size_list[-1],
                                                     include_somatic=1,
                                                     include_germline=2)
                for pep in A_peptides:
                    neoepitopes[pep].append(A_peptides[pep])
                for pep in B_peptides:
                    neoepitopes[pep].append(B_peptides[pep])
                ## WILL NEED TO CHECK FOR DUPLICATES BETWEEN THE TWO COPIES
                transcriptA.reset(reference=True)
                transcriptB.reset(reference=True)
        ## WRITE NEOEPITOPES TO OUTPUT FILE 
    else:
        sys.exit("".join([args.subparser_name, 
                            " is not a valid software mode"]))

