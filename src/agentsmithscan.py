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
import copy
import os
import random
import re
import collections
from operator import itemgetter
from sortedcontainers import SortedDict
from intervaltree import Interval, IntervalTree
import tempfile
import subprocess
import string
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

def gtf_to_cds(gtf_file, dictdir, pickle_it=True):
    """ References cds_dict to get cds bounds for later Bowtie query

        Keys in the dictionary are transcript IDs, while entries are lists of
            relevant CDS/stop codon data
            Data: [chromosome, start, stop, 1(CDS)/0(Stop codon), 
                    +/- strand, reading frame]
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
    if pickle_it:
        pickle_dict = "".join([dictdir, "/", "transcript_to_CDS.pickle"])
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
            start = cds[1]
            stop = cds[2]
            # Interval coordinates are inclusive of start, exclusive of stop
            if stop > start:
                searchable_tree[chrom][start:stop] = transcript_id
            # else:
                # report an error?
    # Write to pickled dictionary
    if pickle_it:
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
    # Interval coordinates are inclusive of start, exclusive of stop
    cds = cds_tree[chrom].search(start)
    for cd in cds:
        if cd.data not in transcript_ids:
            transcript_ids.add(cd.data)
    return transcript_ids

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
    with open(in_vcf, "r") as f:
        for line in f:
            # Preserve header lines with out change
            if line[0:2] == "##":
                header_lines.append(line.strip("\n"))
            # Adjust column header and variant lines
            else:
                tokens = line.strip("\n").split("\t")
                new_line = "\t".join([tokens[0], tokens[1], tokens[2], 
                                        tokens[3], tokens[4], tokens[5], 
                                        tokens[6], tokens[7], tokens[8], 
                                        tokens[10], tokens[9]])
                other_lines.append(new_line)
    # Write new vcf
    with open(out_vcf, "w") as f:
        for line in header_lines:
            f.write(line + "\n")
        for line in other_lines:
            f.write(line + "\n")

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

class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype """

    def __init__(self, bowtie_reference_index, CDS):
        """ Initializes Transcript object.

            This class assumes edits added to a transcript are properly
            phased, consistent, and nonredundant. Most conspicuously, there
            shouldn't be SNVs or insertions among deleted bases.

            bowtie_reference_index: BowtieIndexReference object for retrieving
                reference genome sequence
            CDS: list of all CDS lines for exactly one transcript from GTF;
                a line can be a list pre-split by '\t' or not yet split
        """
        assert len(CDS) > 0
        self.bowtie_reference_index = bowtie_reference_index
        self.intervals = []
        last_chrom, last_strand = None, None
        for line in CDS:
            if type(line) is str: line = line.strip().split('\t')
            try:
                assert last_chrom == line[0]
            except AssertionError:
                if last_chrom is None: pass
            try:
                assert last_strand == line[6]
            except AssertionError:
                if last_strand is None: pass
            # Use exclusive start, inclusive end 0-based coordinates internally
            self.intervals.extend(
                    [int(line[3]) - 2, int(line[4]) - 1]
                )
            last_chrom, last_strand = line[0], line[6]
        # Store edits to coding sequence only
        self.edits = collections.defaultdict(list)
        self.deletion_intervals = []
        self.chrom = last_chrom
        self.rev_strand = (True if last_strand == '-' else False)
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
        self.last_edits = collections.defaultdict(list)
        self.last_deletion_intervals = []
        # Need to sort to bisect_left properly when editing!
        self.intervals.sort()

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

    def edit(self, seq, pos, mutation_type='V'):
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

            No return value.
        """
        if mutation_type == 'D':
            try:
                deletion_size = int(seq)
            except ValueError:
                deletion_size = len(seq)
            self.deletion_intervals.append((pos - 2, pos + deletion_size - 2))
        else:
            self.edits[pos - 1].append((seq, mutation_type))

    def expressed_edits(self, start=None, end=None, genome=True):
        """ Gets expressed set of edits and transcript intervals.

            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified

            Return value: tuple (defaultdict
                                 mapping edits to lists of (seq, mutation_type)
                                 tuples, interval list; this takes
                                 the form of a flattened 
                                 self.deletion_intervals and includes only
                                 those intervals of the transcript that
                                 are expressed)
        """
        if not genome:
            raise NotImplementedError(
                'Retrieving sequence with transcript coordinates not '
                'yet fully supported.'
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
                return ''
        end_index = bisect.bisect_left(self.intervals, end)
        if not (end_index % 2):
            # end should be end of CDS
            end = self.intervals[end_index - 1]
            end_index -= 1
        intervals = [start - 1] + self.intervals[start_index:end_index] + [end]
        assert len(intervals) % 2 == 0
        # Include only relevant deletion intervals
        relevant_deletion_intervals = []
        if self.deletion_intervals:
            sorted_deletion_intervals = sorted(self.deletion_intervals)
            deletion_intervals = [sorted_deletion_intervals[0][0],
                                    sorted_deletion_intervals[0][1]]
            for i in xrange(1, len(sorted_deletion_intervals)):
                if (sorted_deletion_intervals[i][0] <= deletion_intervals[-1]):
                    deletion_intervals[-2] = min(deletion_intervals[-2],
                                            sorted_deletion_intervals[i][0])
                    deletion_intervals[-1] = max(deletion_intervals[-1],
                                            sorted_deletion_intervals[i][1])
                else:
                    deletion_intervals.extend(
                            list(sorted_deletion_intervals[i])
                        )
            for i in xrange(0, len(deletion_intervals), 2):
                start_index = bisect.bisect_left(intervals,
                                                    deletion_intervals[i])
                end_index = bisect.bisect_left(intervals,
                                                deletion_intervals[i+1])
                if start_index == end_index:
                    if start_index % 2:
                        # Entirely in a single interval
                        relevant_deletion_intervals.extend(
                                deletion_intervals[i:i+2]
                            )
                    # else deletion is entirely outside CDS within start/end
                else:
                    assert end_index > start_index
                    if start_index % 2:
                        pos = deletion_intervals[i]
                    else:
                        pos = intervals[start_index]
                        start_index += 1
                    # deletion_intervals[i] becomes a new end
                    relevant_deletion_intervals.extend(
                            [pos, intervals[start_index]]
                        )
                    if end_index % 2:
                        end_pos = deletion_intervals[i+1]
                    else:
                        end_index -= 1
                        end_pos = intervals[end_index]
                    relevant_deletion_intervals.extend(
                            [intervals[i] for i in
                             xrange(start_index + 1, end_index)]
                        )
                    relevant_deletion_intervals.append(end_pos)
        intervals = sorted(intervals + relevant_deletion_intervals)
        edits = collections.defaultdict(list)
        for pos in self.edits:
            # Add edit if and only if it's in one of the CDSes
            start_index = bisect.bisect_left(intervals, pos)
            for edit in self.edits[pos]:
                if edit[1] == 'V':
                    if start_index % 2:
                        # Add edit if and only if it lies within CDS boundaries
                        edits[pos].append(edit)
                elif edit[1] == 'I':
                    if start_index % 2 or pos == intervals[start_index]:
                        '''An insertion is valid before or after a block'''
                        edits[pos].append(edit)
        return (edits, intervals)

    def save(self):
        """ Creates save point for edits.

            No return value.
        """
        self.last_edits = copy.copy(self.edits)
        self.last_deletion_intervals = copy.copy(self.deletion_intervals)

    def seq(self, start=None, end=None, genome=True):
        """ Retrieves transcript sequence between start and end coordinates.

            start: start position (1-indexed, inclusive); None means start of
                transcript
            end: end position (1-indexed, inclusive); None means end of
                transcript
            genome: True iff genome coordinates are specified

            Return value: transcript (sub)sequence
        """
        if end < start: return ''
        # Use 0-based coordinates internally
        if start is None:
            start = self.intervals[0] + 2
        if end is None:
            end = self.intervals[-1] + 1
        if genome:
            # Capture only sequence between start and end
            edits, intervals = self.expressed_edits(start, end,
                                                    genome=True)
            '''Check for insertions at beginnings of intervals, and if they're
            present, shift them to ends of previous intervals so they're
            actually added.'''
            new_edits = copy.copy(edits)
            for i in xrange(0, len(intervals), 2):
                if intervals[i] in edits:
                    assert (len(edits[intervals[i]]) == 1
                                and edits[intervals[i]][0][1] == 'I')
                    if i:
                        new_edits[intervals[i-1]] = new_edits[intervals[i]]
                        del new_edits[intervals[i]]
                    else:
                        intervals = [-1, -1] + intervals
                        new_edits[-1] = new_edits[intervals[i]]
                        del new_edits[intervals[i]]
            seqs = []
            for i in xrange(0, len(intervals), 2):
                seqs.append(
                        self.bowtie_reference_index.get_stretch(
                                self.chrom, intervals[i] + 1,
                                intervals[i + 1] -
                                intervals[i]
                            )
                    )
            # Now build sequence in order of increasing edit position
            i = 1
            pos_group, final_seq = [], []
            for pos in (sorted(new_edits.keys()) + [self.intervals[-1] + 1]):
                if pos > intervals[i]:
                    last_index, last_pos = 0, intervals[i-1] + 1
                    for pos_to_add in pos_group:
                        fill = pos_to_add - last_pos
                        final_seq.append(seqs[(i-1)/2][
                                            last_index:last_index + fill
                                        ])
                        # If no edits, snv is reference and no insertion
                        try:
                            snv = seqs[(i-1)/2][last_index + fill]
                        except IndexError:
                            '''Should happen only for insertions at beginning
                            of sequence.'''
                            assert (i - 1) / 2 == 0 and not seqs[0]
                            snv = ''
                        insertion = ''
                        for edit in new_edits[pos_to_add]:
                            if edit[1] == 'V':
                                snv = edit[0]
                            else:
                                assert edit[1] == 'I'
                                insertion = edit[0]
                        final_seq.extend([snv, insertion])
                        last_index += fill + 1
                        last_pos += fill + 1
                    final_seq.append(seqs[(i-1)/2][last_index:])
                    i += 2
                    try:
                        while pos > intervals[i]:
                            final_seq.append(seqs[(i-1)/2])
                            i += 2
                    except IndexError:
                        if i > len(intervals) - 1:
                            # Done enumerating sequence
                            break
                    pos_group = [pos]
                else:
                    pos_group.append(pos)
            if self.rev_strand:
                return ''.join(final_seq[::-1].translate(
                                    _revcomp_translation_table
                                ))
            return ''.join(final_seq)
        raise NotImplementedError(
            'Retrieving sequence with transcript coordinates not '
            'yet fully supported.'
        )

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

def process_haplotypes(hapcut_output, cds_dict, interval_dict, reference_index,
                        VAF_pos, size_list):
    tumor_peptides = []
    VAFs = []
    with open(hapcut_output, "r") as f:
        block_transcripts = {}
        for line in f:
            if line.startswith('BLOCK'):
                # Skip block header lines
                continue
            elif line[0] == "*":
                # Process all transcripts for the block
                for transcript_ID in block_transcripts:
                    VAF_list = []
                    block_transcripts[transcript_ID].sort(key=itemgetter(1))
                    # Create transcript object and make relevant edits
                    transcript = Transcript(reference_index, 
                                            [[str(chrom), 'blah', 'blah', 
                                                str(start), str(end), 
                                                '.', strand] for (chrom, start, 
                                                                    end, 
                                                                    entry_type, 
                                                                    strand, 
                                                                    frame) in 
                                                cds_dict[transcript_ID]]
                                            )
                    for mutation in block_transcripts[transcript_ID]:
                        if VAF_pos is not None:
                            if "%" in mutation[7].split(":")[VAF_pos]:
                                VAF_list.append(
                                                float(
                                                    mutation[7].split(
                                                                        ":"
                                                        )[VAF_pos].strip(
                                                                        "%"
                                                                        )
                                                    )
                                                )
                        if len(mutation[5]) == len(mutation[6]):
                            mutation_type = 'V'
                        elif len(mutation[5]) < len(mutation[6]):
                            mutation_type = 'I'
                        elif len(mutation[5]) > len(mutation[6]):
                            mutation_type = 'D'
                        else:
                            mutation_type = '?'
                        make_edit = transcript.edit(mutation[3], mutation[1], 
                                        mutation_type=mutation_type)
                        if not make_edit:
                            block_transcripts[transcript_ID].remove(mutation)
                    # Calculate VAF for the transcript
                    if len(VAF_list) == 0:
                        transcript_VAF = "NA"
                    elif len(VAF_list) == 1:
                        transcript_VAF = VAF_list[0]
                    else:
                        VAF_list.sort()
                        if args.VAF_freq_calc == "min":
                            transcript_VAF = VAF_list[0]
                        elif args.VAF_freq_calc == "max":
                            transcript_VAF = VAF_list[-1]
                        elif args.VAF_freq_calc == "mean":
                            transcript_VAF = float(sum(VAF_list))/float(
                                                                len(VAF_list))
                        elif args.VAF_freq_calc == "median":
                            transcript_VAF = median(VAF_list)
                        else:
                            sys.exit(" ".join([args.VAF_freq_calc,
                             " is not a valid VAF calculation option"])
                            )
                    # Get neoepitope sequences for each mutation
                    for mutation in block_transcripts[transcript_ID]:
                        nucleotide_sequence = transcript.seq(
                                    start= mutation[1] - 3*size_list[-1] - 1, 
                                    end=mutation[1] + 3*size_list[-1] - 1, 
                                    genome=True)
                        peptides = neoepitopes(nucleotide_sequence,
                                        reverse_strand=transcript.rev_strand, 
                                        min_size=size_list[0], 
                                        max_size=size_list[-1])
                        for pep in peptides:
                            tumor_peptides.append(pep)
                        for i in range(0, len(peptide_lists[0])):
                            VAFs.append(transcript_VAF)
                # Reset transcript dictionary
                block_transcripts = {}
            else:
                # Add mutation to transcript dictionary for the block
                tokens = line.strip("\n").split("\t")
                overlapping_transcripts = get_transcripts_from_tree(
                                                            tokens[3], 
                                                            tokens[4], 
                                                            interval_dict)
                # For each overlapping transcript, add mutation entry
                # Contains chromosome, position, reference, alternate, allele A,
                #   allele B, genotype line from VCF
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
    return tumor_peptides, VAFs

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
                self.Ycds = gtf_to_cds(self.gtf, "NA", pickle_it=False)
                self.Ytree = cds_to_tree(self.Ycds, "NA", pickle_it=False)
            def test_transcript_to_CDS(self):
                """Fails if dictionary was built incorrectly"""
                self.assertEqual(len(self.Ycds.keys()), 220)
            def test_CDS_tree(self):
                """Fails if dictionary was built incorrectly"""
                self.assertEqual(len(self.Ytree.keys()), 1)
                self.assertEqual(len(self.Ytree["Y"]), 2138)
            def test_transcript_extraction(self):
                """Fails if incorrect transcripts are pulled"""
                self.assertEqual(len(get_transcripts_from_tree(
                                                                "Y", 150860, 
                                                                self.Ytree)), 
                                                                10
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
                combinevcf(self.varscan, self.germline, self.outvcf)
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
        class TestTranscript(unittest.TestCase):
            """Tests transcript object construction"""
            def setUp(self):
                """Sets up files for testing"""
                self.fasta = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'ref.fasta'
                                )
                self.ref_prefix = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'ref'
                                )
                with open(self.fasta, "w") as f:
                    f.write(">1 dna:chromosome\n")
                    f.write("ACGCCCGTGACTTATTCGTGTGCAGACTAC\n")
                    f.write("ATGCCCGTGCCGAATTCGTGTCCCCGCTAC\n")
                    f.write("AATGCCCGTGCCGATTTGAAACCCCGCTAC\n")
                subprocess.call(["bowtie-build", self.fasta, self.ref_prefix])
                self.reference_index = bowtie_index.BowtieIndexReference(
                                                                self.ref_prefix)
                self.CDS = ["1", 'blah', 'blah', "31", "75", '.', "+"]
                self.stop = ["1", 'blah', 'blah', "76", "78", '.', "+"]
                self.transcript = Transcript(self.reference_index, 
                                                [self.CDS, self.stop])
            def test_irrelevant_edit(self):
                """Fails if edit is made for non-CDS/stop position"""
                self.transcript.edit("G", 1)
                relevant_edits = self.transcript.expressed_edits()
                self.assertEqual(self.transcript.edits[0], [("G", "V")])
                self.assertEqual(relevant_edits[0], {})
                self.assertEqual(relevant_edits[1], [29, 74, 74, 77])
            def test_relevant_edit(self):
                """Fails if edit is not made for CDS position"""
                self.transcript.edit("G", 34)
                relevant_edits = self.transcript.expressed_edits()
                self.assertEqual(self.transcript.edits[33], [("G", "V")])
                self.assertEqual(relevant_edits[0][33], [("G", "V")])
                self.assertEqual(relevant_edits[1], [29, 74, 74, 77])
            def test_reset_to_reference(self):
                """Fails if transcript is not reset to reference"""
                self.transcript.edit("G", 34)
                self.transcript.reset(reference=True)
                self.assertEqual(self.transcript.edits, {})
            def test_edit_and_save(self):
                """Fails if edits aren't saved"""
                self.transcript.edit("G", 34)
                self.transcript.edit("CCC", 60, mutation_type="D")
                self.transcript.save()
                self.assertEqual(self.transcript.last_edits[33], [("G", "V")])
                self.assertEqual(self.transcript.last_deletion_intervals,
                                    [(59, 62)])
            def test_reset_to_save_point(self):
                """Fails if new edit not erased or old edits not retained"""
                self.transcript.edit("G", 34)
                self.transcript.edit("CCC", 60, mutation_type="D")
                self.transcript.save()
                self.transcript.edit("C", 36)
                self.assertEqual(self.transcript.edits[35], [("C", "V")])
                self.transcript.reset(reference=False)
                self.assertNotIn(35, self.transcript.edits)
                self.assertEqual(self.transcript.last_edits[33], [("G", "V")])
                self.assertEqual(self.transcript.last_deletion_intervals,
                                    [(59, 62)])
                self.assertNotEqual(self.transcript.edits, {})
            def test_SNV_seq(self):
                """Fails if SNV is edited incorrectly"""
                self.transcript.edit("G", 31)
                seq1 = self.transcript.seq()
                seq2 = self.transcript.seq(31, 36)
                self.assertEqual(len(seq1), 48)
                self.assertEqual(len(seq2), 6)
                self.assertEqual(seq1, 
                            "GTGCCCGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG")
                self.assertEqual(seq2, "GTGCCC")
            def test_inside_indel(self):
                """Fails if indel within CDS is inserted incorrectly"""
                self.transcript.edit("Q", 35, mutation_type="I")
                self.assertEqual(self.transcript.edits[34], [("Q", "I")])
                seq1 = self.transcript.seq()
                seq2 = self.transcript.seq(31, 36)
                self.assertEqual(len(seq1), 49)
                self.assertEqual(len(seq2), 7)
                self.assertEqual(seq1, 
                            "ATGCCQCGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG")
                self.assertEqual(seq2, "ATGCCQC")
            def test_adjacent_indel(self):
                """Fails if indel right before CDS is inserted incorrectly"""
                self.transcript.edit("Q", 30, mutation_type="I")
                self.assertEqual(self.transcript.edits[29], [("Q", "I")])
                seq1 = self.transcript.seq()
                seq2 = self.transcript.seq(31, 36)
                self.assertEqual(len(seq1), 49)
                #self.assertEqual(len(seq2), 6)
                self.assertEqual(seq1, 
                            "QATGCCCGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG")
                #self.assertEqual(seq2, "QATGCC")
                ### SEQ2 RETURNS ATGCCCQATGCCC
            def test_deletion(self):
                """Fails if deletion is made incorrectly"""
                self.transcript.edit(5, 34, mutation_type="D")
                self.assertEqual(self.transcript.deletion_intervals, [(33, 38)])
                seq1 = self.transcript.seq()
                seq2 = self.transcript.seq(31, 36)
                self.assertEqual(len(seq1), 43)
                self.assertEqual(len(seq2), 6)
                self.assertEqual(seq1, 
                                "ATGGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG")
                self.assertEqual(seq2, "ATGGCC")
                ### SEQ2 can't be obtained
            def tearDown(self):
                """Removes temporary files"""
                ref_remove = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'ref*ebwt'
                                )
                subprocess.call(["rm", ref_remove])
                subprocess.call(["rm", self.fasta])
        class TestHAPCUT2Processing(unittest.TestCase):
            """Tests proper processing of HAPCUT2 files"""
            def setUp(self):
                """Sets up input files and dictionaries to use for tests"""
                self.gtf = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.gtf'
                                )
                self.Ycds = gtf_to_cds(self.gtf, "NA", pickle_it=False)
                self.Ytree = cds_to_tree(self.Ycds, "NA", pickle_it=False)
                self.Yhapcut = os.path.join(
                                    os.path.dirname(
                                            os.path.dirname(
                                                    os.path.realpath(__file__)
                                                )
                                        ), 'test', 'Ychrom.hap.out'
                                )         
            def test_hap_processing(self):
                '''
                """Fails if file is processed incorrectly"""
                Ynorm, Ytum, YVAF = process_haplotypes(self.Yhapcut, self.Ycds, 
                                                    self.Ytree, None, [8,11])
                self.assertEqual([len(Ynorm), len(Ytum), len(YVAF)], [0,0,0])
                #### WRITE TEST FOR CASE WHERE THERE WILL BE EPITOPES ####
                '''
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
        # Create relevant transcript objects and edit with mutations
        tumor_peptides, VAF_list = process_haplotypes(
                                                    args.merged_hapcut2_output,
                                                    cds_dict, interval_dict,
                                                    VAF_pos, size_list
                                                    )
        ## Get binding affinities for neoepitopes
        tumor_affinities = get_affinity(normal_peptides, args.allele)
        ## Find multi-mapping rate of neoepitopes to proteome??
        # Combine neoepitope data into entries and take unique set
        neoepitope_data = zip(tumor_peptides, 
                                tumor_affinities, 
                                VAF_list)
        unique_neoepitopes = set(neoepitope_data)
        ## Prioritize output based on: affinity, VAF, peptide similarity?
        ## Write to output file
    else:
        sys.exit("".join([args.subparser_name, 
                            " is not a valid software mode"]))

