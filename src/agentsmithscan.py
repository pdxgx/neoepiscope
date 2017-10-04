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

def write_neoepitopes(mutation_positions, normal_seq, mutated_seq,
                        reverse_strand=False, min_size=8, max_size=11,
                        output_stream=sys.stdout):
    """ Prints neoepitopes from normal and mutated seqs.

        mutation_positions: list of mutation positions
        normal_seq: normal nucleotide sequence
        mutated_seq: mutated nucelotide sequence
        reverse_strand: True iff strand is -
        min_size: minimum peptide kmer size to write
        max_size: maximum petide kmer size to write

        No return value.
    """
    for normal_kmer, mutated_kmer in zip(
            kmerize_peptide(
                seq_to_peptide(
                    normal_seq, reverse_strand=reverse_strand),
                min_size=min_size,
                max_size=max_size
            ), kmerize_peptide(
                seq_to_peptide(
                    mutated_seq, reverse_strand=reverse_strand),
                min_size=min_size,
                max_size=max_size
            )):
        print >>sys.stdout, (
            '\t'.join([normal_kmer, mutated_kmer, str(mutation_posits)]))

#def get_cds(transcript_id, mutation_pos_list, seq_length_left,
#              seq_length_right, ordered_cds_dict, mutation_dict):
    ''' References cds_dict to get cds Bounds for later Bowtie query.
        transcript_id: (String) Indicates the transcript the mutation
            is located on.
        mutation_pos_list: (int) Mutation's position on chromosome
        seq_length_left: (int) How many bases must be gathered
            to the left of the mutation
        seq_length_right: (int) How many bases must be gathered to
            the right of the mutation
        Return value: List of tuples containing starting indexes and stretch
        lengths within cds boundaries necessary to acquire the complete 
        sequence necessary for peptide kmerization based on the position 
        of a mutation within a chromosome.
    '''
#    return nucleotide_index_list, mutation_dict, bounds_set

class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype """
    def __init__(bowtie_reference_index, CDS):
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

    def reset(reference=False):
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

    def edit(seq, pos, mutation_type='V'):
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

#def find_stop(query_st, trans_id, line_count, cds_dict, chrom, reference_index, mutation_locs, reverse):
    ''' Queries get_cds() and Bowtie to find a stop codon in the case of a phase
            shift (indel)
        query_st: (int)
        trans_id: ()
        line_count: ()
        cds_dict: ()
        chrom: (int) The chromosome that the mutation is located on
        reference_index: ()
        mutation_locs: (Dict) Dictionary mapping sequence position to 
            actual mutation
        reverse: (Bool) True if mutation represented by reverse strand
        Return value:  
    '''


def get_seq(chrom, start, splice_length, reference_index):
    chr_name = "chr" + chrom #proper
    start -= 1 #adjust for 0-based bowtie queries
    try:
        seq = reference_index.get_stretch(chr_name, start, splice_length)
    except KeyError:
        return False
    return seq
    
    
def get_affinity(peptides, allele, method):
	''' Takes in peptides and returns their binding affinities to the specified allele 
			based on some prediction method
		peptides: peptides of interest (list of strings)
		allele: HLA allele to use for binding affinity (string) ### May change to list of equal length to peptides
		method: Program to use for binding affinity (string)
		
		Return value: affinities, a list of binding affinities (strings)
	'''
	### Check if allele is valid
	### Check if method is valid
	#peptides.sort(key=len)
	### Break up list into peptides of same length
	### Use subprocess or similar to call the program for each set
	### Parse output of program and store affinities in a list
	# return affinities
	pass


def go():
    """ Entry point """
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=_help_intro, 
                formatter_class=help_formatter)
    try:
        if args.vcf == '-':
            if sys.stdin.isatty():
                raise RuntimeError('Nothing piped into this script, but input is '
                                   'to be read from stdin')
            else:
                input_stream = sys.stdin
        else:
            input_stream = open(args.vcf, "r")
            line_count = 0
            trans_lines = []
            my_dict = cds_dict
            #print(get_seq("3", 69990, 20, reference_index))
            for line in input_stream:
                line_count += 1
                if not line or line[0] == '#': continue
                vals = line.strip().split('\t')
                info = vals[7]
                tokens = info.strip().split('|')
                (trans_id) = (tokens[6])
                if((len(trans_lines) != 0) and (last_trans != trans_id)):
                    try:
                        direct = orf_dict[last_trans][0][0]
                        cds_list = cds_dict[last_trans]
                    except:
                        print "orf_dict failure"
                        last_trans = trans_id
                        trans_lines = []
                        continue
                    #print "here"
                    if(direct == "-"):
                        trans_lines = list(reversed(trans_lines))
                    begin_line = line_count - len(trans_lines) - 1
                    print "Before kmerize trans", last_trans, len(trans_lines)
                    kmerize_trans(trans_lines, begin_line, last_trans, cds_list, direct)
                    trans_lines = []
                last_trans = trans_id
                trans_lines.append(line)
            try:
                direct = orf_dict[last_trans][0][0]
            except KeyError:
                pass
            begin_line = line_count - len(trans_lines) - 1
            if(direct == "-"):
                trans_lines = list(reversed(trans_lines))
            try:
                cds_list = cds_dict[trans_id]
                kmerize_trans(trans_lines, begin_line, last_trans, cds_list, direct)
            except KeyError:
                pass

    finally:
        if args.vcf != '-':
            input_stream.close()
    try:
        if "," in args.kmer_size:
            (_size_min, _size_max) = args.kmer_size.split(",")
            _size_min = int(_size_min)
            _size_max = int(_size_max)
        else:
            _size_min = int(args.kmer_size)
            _size_max = _size_min
        if (_size_min < 1 or _size_max < 1):
            except ValueError:
                print "Kmer size(s) must be >= 1"
                pass
        if (_size_max < _size_min):
            except ValueError:
                print "Max kmer size cannot be less than min kmer size"
                pass
    except:
        print "Unable to import kmer size from command line parameter, defaulting to 8-11aa kmers"
        _size_min = 8
        _size_max = 11


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', type=str, required=True,
            default='-',
            help='input vcf or "-" for stdin'
        )
    parser.add_argument('-x', '--bowtie-index', type=str, required=True,
            help='path to Bowtie index basename'
        )
    parser.add_argument('-d', '--dicts', type=str, required=True,
            help='input path to pickled dictionaries'
        )
    parser.add_argument('-b', '--bam', action='store_true', required=False,
            default = False, help='T/F bam is used'
        )
    parser.add_argument('-k', '--kmer-size', type=str, required=False,
            default='8,11', help='kmer size for epitope calculation'
        )
    args = parser.parse_args()
    reference_index = bowtie_index.BowtieIndexReference(args.bowtie_index)
    with open(args.dicts, 'rb') as dict_stream:
        (cds_dict, orf_dict,
            exon_dict, exon_orf_dict) = pickle.load(dict_stream)
