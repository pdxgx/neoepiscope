#!/usr/bin/env python
"""
neoscan

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
import defaultdict
import copy
import os
import random
import re
import collections
from operator import itemgetter
from sortedcontainers import SortedDict
from intervaltree import Interval, IntervalTree
import tempfile
#import Hapcut2interpreter as hap

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

def gtf_to_cds(gtf_file, dictdir):
    """ References cds_dict to get cds bounds for later Bowtie query

        Keys in the dictionary are transcript IDs, while entries are lists of
            relevant CDS/stop codon data
            Data: [chromosome, start, stop, 1(CDS)/0(Stop codon), 
                    reading frame, +/1 strand]
        Writes cds_dict as a pickled dictionary

        gtf_file: input gtf file to process
        dictdir: path to directory to store pickled dicts

        Return value: dictionary
    """
    cds_dict = {}
    # Parse GTF to obtain CDS/stop codon info
    with open(gtf_file, "r") as f:
        for line in f:
            if line[0] != '#':
                tokens = line.strip().split('\t')
                if tokens[2] == "CDS" or tokens[2] == "stop_codon":
                    transcript_id = re.sub(
                                r'.*transcript_id \"([A-Z0-9._]+)\"[;].*', 
                                r'\1', tokens[8]
                                )
                    # Create new dictionary entry for new transcripts
                    if transcript_id not in cds_dict:
                        if tokens[2] == "CDS":
                            cds_dict[transcript_id] = [[tokens[0].replace(
                                                                    "chr", ""), 
                                                        int(tokens[3]), 
                                                        int(tokens[4]), 1, 
                                                        tokens[6], 
                                                        int(tokens[7])]]
                        else:
                            cds_dict[transcript_id] = [[tokens[0].replace(
                                                                    "chr", ""), 
                                                        int(tokens[3]), 
                                                        int(tokens[4]), 0, 
                                                        tokens[6], 
                                                        int(tokens[7])]]
                    # Append previous entry for old transcripts
                    else:
                        if tokens[2] == "CDS":
                            cds_dict[transcript_id].append([tokens[0].replace(
                                                                    "chr", ""), 
                                                            int(tokens[3]), 
                                                            int(tokens[4]), 1, 
                                                            tokens[6], 
                                                            int(tokens[7])])
                        else:
                            cds_dict[transcript_id].append([tokens[0].replace(
                                                                    "chr", ""), 
                                                            int(tokens[3]), 
                                                            int(tokens[4]), 0, 
                                                            tokens[6], 
                                                            int(tokens[7])])
    # Sort cds_dict coordinates (left -> right) for each transcript                                
    for transcript_id in cds_dict:
            cds_dict[transcript_id].sort(key=lambda x: x[0])
    # Write to pickled dictionary
    pickle_dict = "".join([dictdir, "/", "transcript_to_CDS.pickle"])
    with open(pickle_dict, "wb") as f:
        pickle.dump(cds_dict, f)
    return cds_dict

def cds_to_tree(cds_dict, dictdir):
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
            start = cds[1]
            stop = cds[2]
            # Interval coordinates are inclusive of start, exclusive of stop
            if stop > start:
                searchable_tree[chrom][start:stop] = transcript_id
            # else:
                # report an error?
    # Write to pickled dictionary
    pickle_dict = "".join([dictdir, "/", "intervals_to_transcript.pickle"])
    with open(pickle_dict, "wb") as f:
        pickle.dump(searchable_tree, f)
    return searchable_tree

def get_transcripts_from_tree(chrom, start, cds_tree):
    """ Uses cds tree to btain transcript IDs from genomic coordinates
            
        chrom: (String) Specify chrom to use for transcript search.
        start: (Int) Specify start position to use for transcript search.
        stop: (Int) Specify ending position to use for transcript search
        cds_tree: (Dict) dictionary of IntervalTree() objects containing
            transcript IDs as function of exon coords indexed by chr/contig ID.
            
        Return value: (set) a set of matching unique transcript IDs.
    """
    transcript_ids = set()
    if stop <= start:
        return transcript_ids
    # Interval coordinates are inclusive of start, exclusive of stop
    cds = cds_tree[chrom].search(start)
    for cd in cds:
        if cd.data not in transcript_ids:
            transcript_ids.add(cd.data)
    return transcript_ids

def combinevcf(vcf1, vcf2, outfile="Combined.vcf"):
    """ Combines VCFs
        
        #### WE NEED TO ADJUST THIS ####
        ## Where are header lines going? ##
        ## Position of tumor vs. normal in somatic vcf is variable ##

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
    markgermline = "".join(['''awk '{print $0, "0"}' ''', vcf2, 
                            ".germlinetemp > ", vcf2, ".germline"])
    marktumor    = "".join(['''awk '{print $0, "1"}' ''', vcf2, 
                            ".tumortemp > ", vcf2, ".tumor"])
    os.system(markgermline)
    os.system(marktumor)
    command = "".join(["cat ", vcf2, ".germline ", vcf2, ".tumor > ", 
                        vcf2, ".combine1"])
    os.system(command)
    command2 = "".join(["sort -k1,1 -k2,2n ", vcf2, ".combine1 > ", 
                        vcf2, ".sorted"])
    os.system(command2)
    command3 = "".join(["cat ", vcf2, ".header ", vcf2, ".sorted > ", 
                        vcf2, ".combine2"])
    os.system(command3)
    cut = "".join(["cut -f1,2,3,4,5,6,7,8,9,10 ", vcf2, 
                    ".combine2 > ", outfile])
    os.system(cut)
    for file in [".tumortemp", ".germlinetemp", ".combine1", ".combine2", 
                    ".sorted", ".tumor", ".germline", ".header"]:
        cleanup = "".join(["rm ", vcf2, file])
        os.system(cleanup)

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

def create_haplotype_dictionary(haplooutfile, freqpos):
    """
        Returns dictionary of haplotype information
        
        haplooutfile: name of haplotype output file
        freqpos: location in genotype field where allele frequency is (0 based)

        outputs: dictionary with key of (chrome, pos) and value of 
            (allele on chrome A, allele on chrome B, reference, 
            variant, allele frequency, phasing score)
    """
    hapfile = open(haplooutfile, "r")
    haplotype_dict = collections.OrderedDict()
    for line in hapfile:
        hapline = line.strip().split()
        if hapline[0] != "********" and hapline[0] != "BLOCK:":
            infoline = hapline[7].strip().split(':')
            if freqpos is None:
                haplotype_dict[(int(hapline[3]), 
                    int(hapline[4]))] = (int(hapline[1]), int(hapline[2]), 
                    hapfileapline[5], hapline[6], 
                    None, float(hapline[10]))
            else:
                haplotype_dict[(int(hapline[3]), 
                    int(hapline[4]))] = (int(hapline[1]), int(hapline[2]), 
                    hapfileapline[5], hapline[6], 
                    float(infoline[freqpos].rstrip('%')), float(hapline[10]))
    return SortedDict(haplotype_dict)

class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype """

    def __init__(self, bowtie_reference_index, CDS):
        """ Initializes Transcript object

            bowtie_reference_index: BowtieIndexReference object for retrieving
                reference genome sequence
            CDS: list of all CDS lines for exactly one transcript from GTF
        """
        self.bowtie_reference_index = bowtie_reference_index
        self.intervals = []
        for line in CDS:
            tokens = line.strip().split('\t')
            self.intervals.extend(
                    [(int(tokens[3]), True), (int(tokens[4]) + 1, True)]
                )
        self.edits = defaultdict(list)
        self.reference_intervals = copy.copy(self.intervals)
        intervals_size = len(self.reference_intervals)
        assert intervals_size >= 2 and intervals_size % 2 == 0
        self.chrom = CDS[0][1]
        '''Assume intervals are nonoverlapping! Uncomment following lines to
        check (slower).'''
        # for i in xrange(1, len(self.intervals)):
        #    if self.intervals[i-1] <= self.intervals[i]:
        #        raise RuntimeError(
        #                ('CDS intervals list '
        #                 '"{}" has overlapping intervals.').format(
        #                            self.intervals
        #                        )
        #            )
        # For retrieving save point
        self.last_edits = None
        self.last_intervals = None

    def reset(self, reference=False):
        """ Resets to last save point or reference (i.e., removes all edits).

            reference: if False, tries to reset to last save point, and if that
                doesn't exist, resets to reference. If True, resets to 
                reference.

            No return value.
        """
        if reference or self.last_edits is None:
            self.edits = []
            self.last_intervals = self.reference_intervals
        else:
            self.edits = self.last_edits
            self.intervals = self.last_intervals

    def edit(self, seq, pos, mutation_type='V'):
        """ Adds an edit to the transcript. 

            seq: sequence to add or delete from reference; for deletions, all
                that matters is this sequence has the same length as the 
                sequence to delete
            pos: 0-based coordinate. For insertions, this is the coordinate 
                directly following the inserted sequence. For deletions, this 
                is the coordinate of the first base of the transcript to be
                deleted. Coordinates are always w.r.t. genome.
            mutation_type: V for SNV, I for insertion, D for deletion

            No return value.
        """
        self.intervals.append((pos, False))
        self.edits[pos].append((seq, mutation_type))

    def edit_freq(self, pos, val):
        """modifies allele freq value at location """
        pass

    def get_freq(start=0, end=None, genome=True):
        """ Retrieves allele frequency list between start and end coordinates
        """
        pass

    def save():
        """ Creates save point for edits.

            No return value.
        """
        self.last_edits = copy.copy(self.edits)
        self.last_intervals = copy.copy(self.intervals)

    def seq(start=None, end=None, genome=False):
        """ Retrieves transcript sequence between start and end coordinates.

            start: start position (0-indexed, inclusive); None means start of
                transcript
            end: end position (0-indexed, exclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified

            Return value: transcript (sub)sequence
        """
        assert end is None or end >= start
        self.intervals.sort()
        # Segregate edits by stretch
        on, edits, edit_group = False, [], []
        for i in xrange(len(self.intervals)):
            if self.intervals[i][1]:
                if on:
                    edits.append(edit_group)
                    edit_group = []
                on = not on
            else:
                edit_group.append(i)
        # Get all CDS stretches
        cds_stretches = [bowtie_reference_index.get_stretch(
                                self.chrom,
                                self.reference_intervals[i][0],
                                self.reference_intervals[i+1][0]
                                    - self.reference_intervals[i][0]
                            ) for i in xrange(
                                        len(self.reference_intervals) - 1)]
        if genome:
            if start is None:
                start = self.reference_intervals[0][0]
            if end is None:
                end = self.reference_intervals[-1][0]
            # Mutate all CDSes between start and end
            edited_cds_stretches, break_outer = [], False
            for i, edit_group in enumerate(edits):
                bounds = (self.reference_intervals[i*2][0],
                            self.reference_intervals[i*2+1][0])
                cds_stretch = cds_stretches[i]
                edited_cds_stretch = []
                cds_start = 0
                for j in edit_group:
                    if self.intervals[j][0] < start:
                        continue
                    if self.intervals[j][0] >= end:
                        # no need to continue
                        break_outer = True
                        break
                    else:
                        pos = self.intervals[j][0]
                        for mutation in self.edits[pos]:
                            if mutation[1] == 'V':
                                edited_cds_stretch.append(

                                    )
                                edited_cds_stretches[
                                    mutation[0] - bounds[0]] = mutation[0]
                            elif mutation[1] == 'I':
                                edited_cds_stretches.append(

                                    )
                                edited_cds_stretches[

                                ]
                if break_outer:
                    break
            started, on, on_at_start = False, False, True
            start_index, end_index = 0, len(intervals) - 1
            for i, point in enumerate(self.intervals):
                if point[1]:
                    on = not on
                if point[0] < start:
                    continue
                elif started:
                    if point[0] >= end:
                        end_index = i
                else:
                    started = True
                    on_at_start = on
                    start_index = i
            # Accumulate transcript sequence
            if start < self.intervals[0][0] and self.intervals[0][1]:
                self.intervals = [(start, True)] + self.intervals[1:]
            elif self.intervals[start_index] != (start, True):
                self.intervals = [(start, True)] + self.intervals[start_index:]
            if end > self.intervals[-1][0] and self.intervals[-1][1]:
                self.intervals = self.intervals[:-1] + [(end, True)]
            elif self.intervals[end_index] != (end, True):
                self.intervals = self.intervals[:end_index] + [(end, True)]
            seq = []
            on, i, start = on_at_start, 0, None
            while True:
                if not self.intervals[i][1]:
                    seq.append(self.bowtie_reference_index.get_stretch(
                                                    start, self.intervals[i][]
                                                )
                        )
                if on:
                    reference_index.get_stretch()
                on = not on
            reference_index.get_stretch()
            reference_index.get_stretch(start, self.intervals[i])

        raise NotImplementedError(
            'Retrieving sequence with transcript coordinates not '
            'yet supported.'
        )

def phase_mutations(transcript1, transcript2, hapdict, chromosome, start, end):
    '''
        Applies phased mutations to 2 transcripts
        
        transcript1: transcript object for first chromosome
        transcript2: transcript object for second chromosome
        hapdict: haplotype dictionary
        chromosome: chromosome of transcript
        start: start of applying mutations
        end: end of applying mutations

        Return values: transcript1, transcript2 (edited)
    '''
    startindex = hapdict.bisect((chromosome, start)) - 1
    endindex = hapdict.bisect((chromosome, end)) + 1
    for x in range(startindex, endindex):
        mute_type = 'V'
        if start <= hapdict.items()[x][0][1] <= end:
            # Determine mutation type - default = SNV ("V")
            if len(hapdict.items()[x][1][3]) > len(hapdict.items()[x][1][2]):
                mute_type = 'I'
            elif len(hapdict.items()[x][1][3]) < len(hapdict.items()[x][1][2]):
                mute_type = 'D'
            # Edit transcripts based on mutation type
            if mute_type = 'V':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][3], 
                                        hapdict.items()[x][0][1], mute_type)
                    transcript1.edit_freq(hapdict.items()[x][0][1], 
                                            hapdict.items()[x][1][5])
                else:
                    transcript2.edit(hapdict.items()[x][1][3], 
                                        hapdict.items()[x][0][1], mute_type)
                    transcript2.edit_freq(hapdict.items()[x][0][1], 
                                            hapdict.items()[x][1][5])
            if mute_type = 'I':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][3][1:], 
                                        hapdict.items()[x][0][1]+1, mute_type)
                    transcript1.edit_freq(hapdict.items()[x][0][1]+1, 
                                            hapdict.items()[x][1][5])
                else:
                    transcript2.edit(hapdict.items()[x][1][3][1:], 
                                        hapdict.items()[x][0][1]+1, mute_type)
                    transcript2.edit_freq(hapdict.items()[x][0][1]+1, 
                                            hapdict.items()[x][1][5])
            elif mute_type = 'D':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][2][1:], 
                                        hapdict.items()[x][0][1]+1, mute_type)
                    transcript1.edit_freq(hapdict.items()[x][0][1]+1, 
                                            hapdict.items()[x][1][5])
                else:
                    transcript2.edit(hapdict.items()[x][1][2][1:], 
                                        hapdict.items()[x][0][1]+1, mute_type)
                    transcript2.edit_freq(hapdict.items()[x][0][1]+1, 
                                            hapdict.items()[x][1][5])
            else:
                raise ValueError(" ".join(["Invalid mutation type encountered:",
                                            mute_type]))
    return transcript1, transcript2

def seq_to_peptide(seq, reverse_strand=False):
    """ Translates nucleotide sequence into peptide sequence.

        All codons including and after stop codon are recorded as X's.

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
    for j in xrange(i + 3, seq_size - seq_size % 3, 3):
        peptide.append('X')
    return ''.join(peptide)

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

def neoepitopes(normal_seq, mutated_seq,
                        reverse_strand=False, min_size=8, max_size=11,):
    """ Finds neoepitopes from normal and mutated seqs.

        normal_seq: normal nucleotide sequence
        mutated_seq: mutated nucelotide sequence
        reverse_strand: True iff strand is -
        min_size: minimum peptide kmer size to write
        max_size: maximum petide kmer size to write

        Return value: List of tuples (normal_kmer, mutated_kmer)
    """
    return zip(kmerize_peptide(
        seq_to_peptide(
            normal_seq, reverse_strand=reverse_strand),
        min_size=min_size,
        max_size=max_size
    ), kmerize_peptide(
        seq_to_peptide(
            mutated_seq, reverse_strand=reverse_strand),
        min_size=min_size,
        max_size=max_size))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro, 
                                        formatter_class=help_formatter)
    subparsers = parser.add_subparsers(help=(
                                    'subcommands; add "-h" or "--help" '
                                    'after a subcommand for its parameters'
                                ), dest='subparser_name')
    index_parser = subparsers.add_parser('index',
                                        help=('produces pickled dictionaries '
                                        'linking transcripts to intervals and '
                                        ' CDS lines in a GTF'), 
                                        dest='subparser_name')
    prep_parser = subparsers.add_parser('prep',
                                        help=('combines HAPCUT2 output with '
                                              'unphased variants for call mode'),
                                        dest='subparser_name')
    merge_parser = subparsers.add_parser('merge',
                                        help=('merges germline and somatic '
                                            'VCFS for combined mutation '
                                            'phasing with HAPCUT2'), 
                                        dest='subparser_name')
    call_parser = subparsers.add_parser('call', help='calls neoepitopes', 
                                        dest='subparser_name')
    prep_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input VCF'
        )
    prep_parser.add_argument('-c', '--hapcut2-output', type=str, required=True,
            help='path to output file of HAPCUT2 run on input VCF'
        )
    prep_parser.add_argument('-o', '--output', type=str, required=True,
            help='path to output file to be input to call mode'
        )
    index_parser.add_argument('-g', '--gtf', type=str, required=False,
            help='input path to GTF file'
        )  
    index_parser.add_argument('-d', '--dicts', type=str, required=False,
            help='output path to pickled CDS dictionary'
        )
    merge_parser.add_argument('-g', '--germline', type=str, required=False,
            help='input path to germline VCF'
        )
    merge_parser.add_argument('-s', '--somatic', type=str, required=False,
            help='input path to somatic VCF'
        )
    merge_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to combined VCF'
        )
    call_parser.add_argument('-x', '--bowtie-index', type=str, required=True,
            help='path to Bowtie index basename'
        )
    call_parser.add_argument('-d', '--dicts', type=str, required=False,
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
    args = parser.parse_args()
    
    if args.subparser_name == 'index':
        cds_dict = gtf_to_cds(args.gtf, args.dicts)
        cds_to_tree(cds_dict, args.dicts)
        # FM indexing of proteome??
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
                    if (tokens[3], tokens[4]) not in phased[
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
    elif args.subparser_name == 'merge':
        pass
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
                            tokens = line.split('\t'):
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
        
        # Obtain VAF frequency VCF position and parse hapcut2 output
        VAF_pos = get_VAF_pos(args.vcf)
        hap_dict = create_haplotype_dictionary(args.hapcut2_output, VAF_pos)

        # Obtain peptide sizes for kmerizing peptides
        if ',' in args.kmer_size:
            size_list = args.kmer_size.split(',')
            size_list.sort()
        else:
            size_list = [args.kmer_size]

        # For retrieving genome sequence
        reference_index = bowtie_index.BowtieIndexReference(args.bowtie_index)
        
        # Find transcripts that mutations overlap 
        # Create relevant transcript objects and store in dictionary
        # Store coordinates of unphased mutations
        transcript_dict = {}
        unphased_mutations = {}
        with open(args.merged_hapcut2_output, "r") as f:
            block_mutations = []
            block_transcripts = {}
            for line in f:
                if line.startswith('BLOCK'):
                    continue
                elif line[0] == "*":
                    for transcript_ID in block_transcripts:
                        block_transcripts[transcript_ID].sort(key=itemgetter(1))
                        transcript = Transcript(reference_index, 
                                                cds_dict[transcript_ID])
                        base_sequence = transcript.seq()
                    
                    block_transcripts = {}
                else:
                    tokens = line.strip("\n").split("\t")
                    overlapping_transcripts = get_transcripts_from_tree(
                                                                tokens[3], 
                                                                locus[4], 
                                                                interval_dict)
                    for transcript in overlapping_transcripts:
                        if transcript not in block_transcripts:
                            block_transcripts[transcript] = [[tokens[3], 
                                                                tokens[4], 
                                                                tokens[5], 
                                                                tokens[6], 
                                                                tokens[1], 
                                                                tokens[2], 
                                                                tokens[7]]]
                        else:
                            block_transcripts[transcript].append([tokens[3], 
                                                                tokens[4], 
                                                                tokens[5], 
                                                                tokens[6], 
                                                                tokens[1], 
                                                                tokens[2], 
                                                                tokens[7]])







            for line in f:
                if line[0] != "#":
                    tokens = line.strip("\n").split("\t")
                    locus = (tokens[0].replace("chr", ""), int(tokens[1]))
                    overlapping_transcripts = get_transcripts_from_tree(
                                                                locus[0], 
                                                                locus[1], 
                                                                interval_dict)
                    for transcript in overlapping_transcripts:
                        if transcript not in transcript_dict:
                            transcript_dict[transcript] = Transcript(
                                                            reference_index, 
                                                            cds_dict[transcript]
                                                            )
                        if locus not in hap_dict:
                            if locus not in unphased_mutations:
                                unphased_mutations[locus] = [transcript]
                            else:
                                unphased_mutations[locus].append(transcript)


        ## PHASE MUTATIONS FROM HAPCUT2 ##
        ## MAKE COMBINATORIAL EDITS FOR UNPHASED MUTATIONS ##
        
        ## TRANSLATE PEPTIDES SEQUENCES ##

        ## KMERIZE PEPTIDES ##
        # Call neoepitopes function to produce paired_peptides list
        paired_peptides =[] # unique set of tuples [(norm, tum), (norm, tum)]
        peptide_lists = map(list, zip(*paired_peptides)) # normal then tumor

        ## Get binding affinities for neoepitopes and paired normal epitopes
        normal_affinities = get_affinity(peptides_lists[0], args.allele)
        tumor_affinities = get_affinity(peptides_lists[1], args.allele)


        ## Find multi-mapping rate of neoepitopes to proteome?
        ## Prioritize output based on: affinity, VAF, peptide similarity?

