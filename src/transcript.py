import bowtie_index
import collections
import copy
import bisect
import string
import re
import pickle
from intervaltree import Interval, IntervalTree
from operator import itemgetter

revcomp_translation_table = string.maketrans('ATCG', 'TAGC')

def custom_bisect_left(a, x, lo=0, hi=None, getter=0):
    """ Same as bisect.bisect_left, but compares only index "getter"

        See bisect_left source for more info.
    """

    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid][getter] < x: lo = mid+1
        else: hi = mid
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
    return [item for sublist in
                [[peptide[i:i+size] for i in xrange(peptide_size - size + 1)]
                    for size in xrange(min_size, max_size + 1)]
            for item in sublist if 'X' not in item]

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
            return ''
    seq_size = len(seq)
    peptide = []
    for i in xrange(0, seq_size - seq_size % 3, 3):
        codon = _codon_table[seq[i:i+3]]
        peptide.append(codon)
        if codon == 'X':
            break
    # for j in xrange(i + 3, seq_size - seq_size % 3, 3):
        # peptide.append('X')
    return ''.join(peptide)


class Transcript(object):
    """ Transforms transcript with edits (SNPs, indels) from haplotype. """

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
        self.start_codon, self.stop_codon = None, None
        last_chrom, last_strand = None, None
        for line in CDS:
            if type(line) is str: line = line.strip().split('\t')
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
                else: raise
            # Use exclusive start, inclusive end 0-based coordinates internally
            if line[2] == "exon":
                self.intervals.extend(
                        [int(line[3]) - 2, int(line[4]) - 1]
                    )
            elif line[2] == "start_codon":
                self.start_codon = int(line[3]) - 1
            elif line[2] == "stop_codon":
                self.stop_codon = int(line[3]) - 1
            else:
                raise NotImplementedError(
                                    'GTF sequence type not currently supported'
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
        if self.start_codon:
            self.start_codon_index = bisect.bisect_left(self.intervals, 
                                                        self.start_codon)
        else:
            self.start_codon_index = None
        if self.stop_codon:
            self.stop_codon_index = bisect.bisect_left(self.intervals, 
                                                        self.stop_codon)
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

    def edit(self, seq, pos, mutation_type='V', mutation_class='S', vaf=None):
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
            No return value.
        """
        if mutation_type == 'D':
            try:
                deletion_size = int(seq)
            except ValueError:
                deletion_size = len(seq)
                self.deletion_intervals.append(
                        (pos - 2, pos + deletion_size - 2, mutation_class, 
                            (pos, seq, mutation_type), vaf)
                    )
            else:
                self.deletion_intervals.append(
                        (pos - 2, pos + deletion_size - 2, mutation_class, 
                            (pos, self.bowtie_reference_index.get_stretch(
                                    self.chrom, pos - 1, 
                                    pos + deletion_size + 1  - pos - 1
                                ), mutation_type), vaf)
                    )
        elif mutation_type == 'I': 
            self.edits[pos - 1].append((seq, mutation_type, mutation_class, 
                                        (pos, seq, mutation_type), vaf))
        elif mutation_type == 'V':
            reference_seq = self.bowtie_reference_index.get_stretch(
                                            self.chrom, pos - 1, len(seq))
            self.edits[pos - 1].append((seq, mutation_type, mutation_class, 
                                        (pos, reference_seq, mutation_type), 
                                        vaf))
        else:
            raise NotImplementedError('Mutation type not yet implemented')

    def expressed_edits(self, start=None, end=None, genome=True, 
                                include_somatic=True, include_germline=True):
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
        relevant_deletion_intervals, edits = [], collections.defaultdict(list)
        sorted_deletion_intervals = [
                interval for interval in self.deletion_intervals
                if (interval[2] == 'S' and include_somatic or
                    interval[2] == 'G' and include_germline)
            ]
        if sorted_deletion_intervals:
            sorted_deletion_intervals.sort(key=itemgetter(0, 1))
            deletion_intervals = [(sorted_deletion_intervals[0][0],
                                   sorted_deletion_intervals[0][2], 
                                   sorted_deletion_intervals[0][3],
                                   sorted_deletion_intervals[0][4]),
                                  (sorted_deletion_intervals[0][1],
                                   sorted_deletion_intervals[0][2],
                                   sorted_deletion_intervals[0][3],
                                   sorted_deletion_intervals[0][4])]
            for i in xrange(1, len(sorted_deletion_intervals)):
                if (sorted_deletion_intervals[i][0]
                    <= deletion_intervals[-1][0]):
                    deletion_intervals[-2] = min(deletion_intervals[-2],
                                            (sorted_deletion_intervals[i][0],
                                             sorted_deletion_intervals[i][2],
                                             sorted_deletion_intervals[i][3],
                                             sorted_deletion_intervals[i][4]),
                                            key=itemgetter(0))
                    deletion_intervals[-1] = max(deletion_intervals[-1],
                                            (sorted_deletion_intervals[i][1],
                                             sorted_deletion_intervals[i][2],
                                             sorted_deletion_intervals[i][3],
                                             sorted_deletion_intervals[i][4]),
                                            key=itemgetter(0))
                else:
                    deletion_intervals.extend(
                            [(sorted_deletion_intervals[i][0],
                                sorted_deletion_intervals[i][2], 
                                sorted_deletion_intervals[i][3],
                                sorted_deletion_intervals[i][4]),
                             (sorted_deletion_intervals[i][1],
                                sorted_deletion_intervals[i][2],
                                sorted_deletion_intervals[i][3],
                                sorted_deletion_intervals[i][4])]
                        )
            for i in xrange(0, len(deletion_intervals), 2):
                start_index = bisect.bisect_left(intervals,
                                                    deletion_intervals[i][0])
                end_index = bisect.bisect_left(intervals,
                                                deletion_intervals[i+1][0])
                if start_index == end_index:
                    if start_index % 2:
                        # Entirely in a single interval
                        relevant_deletion_intervals.extend(
                                deletion_intervals[i:i+2]
                            )
                    # else deletion is entirely outside CDS within start/end
                elif (start_index == (end_index - 1) and 
                        deletion_intervals[i][0] == intervals[start_index]):
                    relevant_deletion_intervals.extend(
                                deletion_intervals[i:i+2]
                            )
                else:
                    assert end_index > start_index
                    if (start_index % 2 or 
                        deletion_intervals[i][0] == intervals[start_index]):
                        pos = deletion_intervals[i]
                    else:
                        pos = (intervals[start_index], 'R', tuple(), None)
                        start_index += 1
                    # deletion_intervals[i] becomes a new end
                    relevant_deletion_intervals.extend(
                            [pos, (intervals[start_index], 'R', tuple(), None)]
                        )
                    if end_index % 2:
                        end_pos = deletion_intervals[i+1]
                        relevant_deletion_intervals.extend(
                            [(intervals[i], 'R', tuple(), None) for i in
                             xrange(start_index + 1, end_index)]
                        )
                    else:
                        end_pos = (intervals[end_index - 1], 'R', tuple(), None)
                        relevant_deletion_intervals.extend(
                                [(intervals[i], 'R', tuple(), None) for i in
                                 xrange(start_index, end_index)]
                            )
                    relevant_deletion_intervals.append(end_pos)
        intervals = sorted([(interval, 'R', tuple(), None) for interval in 
                            intervals] + relevant_deletion_intervals)
        edits = collections.defaultdict(list)
        for pos in self.edits:
            # Add edit if and only if it's in one of the CDSes
            start_index = custom_bisect_left(intervals, pos)
            for edit in self.edits[pos]:
                if (include_somatic and edit[2] == 'S'
                        or include_germline and edit[2] == 'G'):
                    if edit[1] == 'V':
                        if start_index % 2 and edit[3][1] != edit[0]:
                            # Add edit if and only if it lies within bounds
                            edits[pos].append(edit)
                    elif edit[1] == 'I':
                        if start_index % 2 or pos == intervals[start_index][0]:
                            # An insertion is valid before or after a block
                            edits[pos].append(edit)
        # Remove empty intervals
        intervals = [intervals[i] for i in xrange(len(intervals))
                         if (i % 2
                             and intervals[i][0] != intervals[i-1][0]
                             or i % 2 == 0
                             and intervals[i+1][0] != intervals[i][0])]
        # Only associate one end of a deletion interval with deletion
        #   to prevent including it multiple times
        adjusted_intervals = [intervals[0]]
        deletion_data = []
        if intervals[0][1] != 'R':
            deletion_data.append(intervals[0][2])
        for i in xrange(1, len(intervals)):
            if intervals[i][1] == 'R':
                adjusted_intervals.append(intervals[i])
            else:
                if intervals[i][2] not in deletion_data:
                    adjusted_intervals.append(intervals[i])
                    deletion_data.append(intervals[i][2])
                else:
                    adjusted_intervals.append((intervals[i][0], 'R', 
                                            tuple(), None))
                # Adjust mutation class to reflect hybrid mutation if needed
                if (intervals[i][1] != intervals[i-1][1] and i % 2 and 
                        intervals[i][1] != 'R' and intervals[i-1][1] != 'R'):
                    mutation_class = ''.join([intervals[i-1][1], 
                                                intervals[i][1]])
                else:
                    mutation_class = intervals[i-1][1]
                # Adjust mutation data to reflect hybrid mutation if needed
                if (intervals[i-1][2] != intervals[i][2] and 
                                intervals[i][2] not in deletion_data and
                                intervals[i-1] != 'R'):
                    mutation_data = []
                    mutation_data.append(intervals[i-1][2])
                    mutation_data.append(intervals[i][2])
                else:
                    mutation_data = intervals[i-1][2]
                # Update previous interval
                adjusted_intervals[i-1] = (adjusted_intervals[i-1][0],
                                                mutation_class,
                                                mutation_data,
                                                adjusted_intervals[i-1][3])
        print intervals
        print
        print adjusted_intervals
        print
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
        if (not (pos_index % 2) or not self.start_codon_index or 
            not self.stop_codon_index):
            # We're outside exon sequence
            return None
        if self.rev_strand:
            if pos_index == self.start_codon_index:
                # Within the same interval as the start codon
                if pos > self.start_codon + 2:
                    # Outside coding sequence
                    return None
                return ((self.start_codon + 2 - pos) % 3)
            else:
                if pos > self.start_codon or pos < self.stop_codon:
                    return None
                seq_length = ((self.intervals[pos_index] - pos + 1) + 
                              (self.start_codon + 2 - 
                                self.intervals[self.start_codon_index - 1]) + 
                              sum([self.intervals[i+1] - self.intervals[i]
                                for i in xrange(pos_index + 1, 
                                                self.start_codon_index - 1, 
                                                2)]))
                return (seq_length - 1) % 3
        else:
            if pos_index == self.start_codon_index:
                if pos < self.start_codon:
                    return None
                return (pos - self.start_codon) % 3
            else:
                if pos < self.start_codon or pos > self.stop_codon:
                    return None
                seq_length = ((pos - self.intervals[pos_index - 1]) + 
                              (self.intervals[self.start_codon_index] - 
                                self.start_codon + 1) + 
                              sum([self.intervals[i+1] - self.intervals[i]
                                for i in xrange(self.start_codon_index + 1,
                                        pos_index - 1, 2)]))
                return ((seq_length - 1) % 3)

    def _seq_append(self, seq_list, seq, mutation_class,
                        mutation_info, vaf, position):
        """ Appends mutation to seq_list, merging successive mutations.
            seq_list: list of tuples (sequence, type) where type is one
                of R, G, or S (for respectively reference, germline edit, or
                somatic edit). Empty sequence means there was a deletion.
            seq: seq to add
            mutation_class: S for somatic, G for germline, R for reference
            mutation_info: tuple containing (1 based mutation position from 
                vcf, mutation sequence, and mutation type)
            vaf: variant allele frequency
            position: 1-based genomic position of first base added
            No return value; seq_list is merely updated.
        """
        try:
            condition = seq_list[-1][1] == mutation_class
        except IndexError:
            # Add first item in seq_list
            assert not seq_list
            if seq or mutation_class != 'R':
                if isinstance(mutation_info, list):
                    seq_list.append((seq, mutation_class,
                                 [mutation_info[i] for i in xrange(0, 
                                                        len(mutation_info))], 
                                 [] if vaf is None else [vaf], position))
                else:
                    seq_list.append((seq, mutation_class,
                                 [mutation_info], 
                                 [] if vaf is None else [vaf], position))
            return
        if condition:
            if self.rev_strand:
                adjusted_position = position
            else:
                adjusted_position = seq_list[-1][4]
            if isinstance(mutation_info, list):
                adjusted_mutation_info = (seq_list[-1][2] + 
                                        [x for x in mutation_info 
                                        if x not in seq_list[-1][2]])
            elif mutation_info not in seq_list[-1][2]:
                adjusted_mutation_info = seq_list[-1][2] + [mutation_info]
            else:
                adjusted_mutation_info = seq_list[-1][2]
            seq_list[-1] = (seq_list[-1][0] + seq, mutation_class,
                            adjusted_mutation_info, 
                            seq_list[-1][3]
                            if vaf is None or vaf in seq_list[-1][3]
                            else (seq_list[-1][3] + [vaf]), adjusted_position)
        elif seq or mutation_class != 'R':
            if isinstance(mutation_info, list):
                seq_list.append((seq, mutation_class,
                                 [mutation_info[i] for i in xrange(0, 
                                                        len(mutation_info))], 
                                 [] if vaf is None else [vaf], position))
            else:
                seq_list.append((seq, mutation_class,
                                 [mutation_info], 
                                 [] if vaf is None else [vaf], position))

    def annotated_seq(self, start=None, end=None, genome=True, 
                                include_somatic=True, include_germline=True):
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
                mutation information, variant allele frequency, position),
                where sequence is a segment of sequence of the (possibly)
                mutated transcript, mutation class is one of {'G', 'S', 'R'},
                where 'G' denotes germline, 'S' denotes somatic, and 'R'
                denotes reference sequence, mutation information is the
                tuple (1-based position of {first base of deletion,
                base before insertion, SNV},
                {deleted sequence, inserted sequence, reference base},
                {'D', 'I', 'V'}) , and position is the 1-based position
                of the first base of sequence.
        """
        if end < start: return ''
        # Use 0-based coordinates internally
        if start is None:
            start = self.intervals[0] + 2
        if end is None:
            end = self.intervals[-1] + 1
        if genome:
            # Capture only sequence between start and end
            edits, intervals = self.expressed_edits(start, end, genome=True, 
                                            include_somatic=include_somatic, 
                                            include_germline=include_germline)
            '''Check for insertions at beginnings of intervals, and if they're
            present, shift them to ends of previous intervals so they're
            actually added.'''
            new_edits = copy.copy(edits)
            i = 0
            while i < len(intervals):
                if intervals[i][0] in edits and i:
                    assert (len(edits[intervals[i][0]]) == 1
                                and edits[intervals[i][0]][0][1] == 'I')
                    new_edits[
                        intervals[i-1][0]] = new_edits[intervals[i][0]]
                    del new_edits[intervals[i][0]]
                    '''Code below would add insertion to first block,
                    but no more
                    if i do the above,
                    else:
                        intervals = [(-1, 'R', [], []), 
                                     (-1, 'R', [], [])] + intervals
                        # Have to add 2 because we modified intervals above
                        i += 2
                        new_edits[-1] = new_edits[intervals[i][0]]
                        del new_edits[intervals[i][0]]'''
                i += 2
            seqs = []
            for i in xrange(0, len(intervals), 2):
                seqs.append(
                        (self.bowtie_reference_index.get_stretch(
                                self.chrom, intervals[i][0] + 1,
                                intervals[i + 1][0] -
                                intervals[i][0]), 
                        (intervals[i][0] + 2, intervals[i+1][0] + 1)
                            )
                    )
            # Now build sequence in order of increasing edit position
            i = 1
            pos_group, final_seq = [], []
            for pos in (sorted(new_edits.keys()) + [self.intervals[-1] + 1]):
                if pos > intervals[i][0]:
                    last_index, last_pos = 0, intervals[i-1][0] + 1
                    for pos_to_add in pos_group:
                        fill = pos_to_add - last_pos
                        if intervals[i-1][1] != 'R':
                            if isinstance(intervals[i-1][2], list):
                                genomic_position = min([x[0] for x 
                                                        in intervals[i-1][2]])
                            else:
                                genomic_position = intervals[i-1][2][0]
                            self._seq_append(final_seq, '', intervals[i-1][1],
                                             intervals[i-1][2],
                                             intervals[i-1][3],
                                             genomic_position)
                        if self.rev_strand:
                            self._seq_append(final_seq, seqs[(i-1)/2][0][
                                            last_index:last_index + fill
                                        ], 'R', tuple(), None,
                                        seqs[(i-1)/2][1][0] + last_index + fill - 1)
                        else:
                            self._seq_append(final_seq, seqs[(i-1)/2][0][
                                            last_index:last_index + fill
                                        ], 'R', tuple(),
                                        None, seqs[(i-1)/2][1][0])
                        # If no edits, snv is reference and no insertion
                        try:
                            snv = (seqs[(i-1)/2][0][last_index + fill], 'R', 
                                    tuple(), None, seqs[(i-1)/2][1][0] + fill)
                            print snv
                        except IndexError:
                            '''Should happen only for insertions at beginning
                            of sequence.'''
                            assert (i - 1) / 2 == 0 and not seqs[0][0]
                            snv = ('', 'R', tuple(), None, seqs[(i-1)/2][1][0] + fill)
                        insertion = ('', 'R', tuple(), None, seqs[(i-1)/2][1][0] + fill)
                        for edit in new_edits[pos_to_add]:
                            if edit[1] == 'V':
                                snv = (edit[0], edit[2], edit[3], edit[4], 
                                        edit[3][0])
                            else:
                                assert edit[1] == 'I'
                                insertion = (edit[0], edit[2], edit[3], edit[4], 
                                                edit[3][0])
                        self._seq_append(final_seq, *snv)
                        self._seq_append(final_seq, *insertion)
                        last_index += fill + 1
                        last_pos += fill + 1
                    if intervals[i-1][1] != 'R':
                        if isinstance(intervals[i-1][2], list):
                            genomic_position = min([x[0] for x 
                                                        in intervals[i-1][2]])
                        else:
                            genomic_position = intervals[i-1][2][0]
                        self._seq_append(final_seq, '', intervals[i-1][1], 
                                         intervals[i-1][2], intervals[i-1][3],
                                         genomic_position)
                    ref_to_add = seqs[(i-1)/2][0][last_index:]
                    if ref_to_add:
                        if self.rev_strand:
                            self._seq_append(
                                final_seq, ref_to_add, 'R', tuple(), None, 
                                seqs[(i-1)/2][1][1]
                            )
                        else:
                            self._seq_append(
                                final_seq, ref_to_add, 'R', tuple(), None, 
                                seqs[(i-1)/2][1][0] + last_index
                            )
                    if intervals[i][1] != 'R':
                        if isinstance(intervals[i][2], list):
                            genomic_position = min([x[0] for x 
                                                        in intervals[i][2]])
                        else:
                            genomic_position = intervals[i][2][0]
                        self._seq_append(final_seq, '', intervals[i][1], 
                                            intervals[i][2], intervals[i][3],
                                            genomic_position)
                    i += 2
                    try:
                        while pos > intervals[i][0]:
                            if intervals[i-1][1] != 'R':
                                if isinstance(intervals[i-1][2], list):
                                    genomic_position = min([x[0] for x 
                                                        in intervals[i-1][2]])
                                else:
                                    genomic_position = intervals[i-1][2][0]
                                self._seq_append(
                                        final_seq, '', intervals[i-1][1], 
                                        intervals[i-1][2], intervals[i-1][3],
                                        genomic_position
                                    )
                            if self.rev_strand:
                                self._seq_append(final_seq, seqs[(i-1)/2][0], 
                                                    'R', tuple(), None, 
                                                    seqs[(i-1)/2][1][1])
                            else:
                                self._seq_append(final_seq, seqs[(i-1)/2][0], 
                                                    'R', tuple(), None, 
                                                    seqs[(i-1)/2][1][0])
                            if intervals[i][1] != 'R':
                                if isinstance(intervals[i][2], list):
                                    genomic_position = min([x[0] for x 
                                                        in intervals[i][2]])
                                else:
                                    genomic_position = intervals[i][2][0]
                                self._seq_append(
                                        final_seq, '', intervals[i][1], 
                                        intervals[i][2], intervals[i][3],
                                        genomic_position
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
                return ([(seq[::-1].translate(revcomp_translation_table),
                            mutation_class, orig_seq, vaf, position)
                            for seq, mutation_class, orig_seq, vaf, position in 
                            final_seq][::-1])
            return final_seq
        raise NotImplementedError(
            'Retrieving sequence with transcript coordinates not '
            'yet fully supported.'
        )

    def neoepitopes(self, min_size=8, max_size=11, include_somatic=1, 
        include_germline=2):
        """ Retrieves list of predicted peptide fragments from transcript that 
            include one or more variants.
            min_size: minimum subpeptide length (specified as # of amino acids)
            max_size: maximum subpeptide length (specified as # of amino acids)
            include_somatic: 0 = do not include somatic mutations, 1 = exlude
            somatic mutations from reference comparison, 2 = include somatic
            mutations in both annotated sequence and reference comparison
            include_germline: 0 = do not include germline mutations, 1 = exlude
            germline mutations from reference comparison, 2 = include germline
            mutations in both annotated sequence and reference comparison
            Return value: list of peptides of desired length.
        """
        if include_somatic == include_germline and include_somatic != 1:
            return []
        if max_size < min_size:
            max_size = min_size
        if min_size < 2:
            return []
        def flatten_annotated_seq(annotated_seq):
            sequence = '' # hold flattened nucleotide sequence
            # extract nucleotide sequence from annotated_seq
            for seq in annotated_seq:
                if seq[1] != 'D':
                    sequence += seq[0]
            return sequence
        annotated_seq = self.annotated_seq(include_somatic=include_somatic != 0, 
            include_germline=include_germline != 0)
        sequence = flatten_annotated_seq(annotated_seq)
        # locate position of start codon (first ATG in sequence)
        variant_peps = kmerize_peptide(seq_to_peptide(sequence, 
            reverse_strand=False, require_ATG=True), 
            min_size=min_size, max_size=max_size)
        if variant_peps == []:
            return []
        reference_seq = self.annotated_seq(include_somatic=include_somatic > 1, 
            include_germline=include_germline > 1)
        ref_sequence = flatten_annotated_seq(reference_seq)
        reference_peps = kmerize_peptide(seq_to_peptide(ref_sequence,
            reverse_strand=False, require_ATG=True), 
            min_size=min_size, max_size=max_size)
        neoepitopes = list(set(variant_peps).difference(reference_peps))
        edits = [seq for seq in annotated_seq if seq[1] != 'R']
        # in order to link each variants to a neoepitope, pick off one at a time
        # and test if any neoepitopes are lost in a list diff -- if any lost,
        # then those variants are REQUIRED for the presence of that peptide!!!
        # easy peasy, independent of complexity BUT will be computationally
        # slower (more robust = good though!!)
        # going to require re-engineering annotated_seq() to select individual
        # edits


        # return list of unique neoepitope sequences
        return neoepitopes

    def peptides(self, min_size=8, max_size=11, include_somatic=1, 
        include_germline=2):
        """ Retrieves list of predicted peptide fragments from transcript that 
            include one or more variants.
            min_size: minimum subpeptide length (specified as # of amino acids)
            max_size: maximum subpeptide length (specified as # of amino acids)
            include_somatic: 0 = do not include somatic mutations, 1 = exlude
            somatic mutations from reference comparison, 2 = include somatic
            mutations in both annotated sequence and reference comparison
            include_germline: 0 = do not include germline mutations, 1 = exlude
            germline mutations from reference comparison, 2 = include germline
            mutations in both annotated sequence and reference comparison
            Return value: list of peptides of desired length.
        """
        if include_somatic == include_germline and include_somatic != 1:
            return []
        if max_size < min_size:
            max_size = min_size
        if min_size < 2:
            return []
        annotated_seq = self.annotated_seq(include_somatic=include_somatic, 
            include_germline=include_germline)
        # if no edits to return, then skip all next steps and return []
        # extract nucleotide sequence from annotated_seq
        sequence = ref_sequence = '' # hold flattened nucleotide sequence
#        reference_seq = self.annotated_seq(include_somatic=include_somatic > 1, 
#            include_germline=include_germline > 1)
#        # extract nucleotide sequence from reference_seq
#        for seq in reference_seq:
#            if seq[1] == 'R' or seq[2][0][2] != 'D':
#                ref_sequence += seq[0]
        #
        # some necessary preprocessing here to find the location of the start
        # codon to use -- for majority of cases, this will simply be the
        # original start_codon for reference sequence!  HOWEVER, it is possible
        # that alteration in start codon and/or creation of new one upstream
        # that could influence actual start site for annotated seq . . . for now
        # ignoring all of this and assuming start site is the same
        #
        start = self.start_codon
        ATGs = []
        ATG_counter1 = ATG_counter2 = 0
        ATG_limit = 2
        strand = 1 - self.rev_strand * 2
        coding_start = ref_start = -1
        counter = ref_counter = 0 # hold edited transcript level coordinates
        seq_previous = []
        # the ATG stuff here is working retrospectively, so may need to run loop
        # one additional time at end to ensure capture start (i.e. what if it occurs
        # during only or final segment in annotated_seq???)
        for seq in annotated_seq:
            # build pairwise list of 'ATG's from annotated_seq and reference
            ATG1 = sequence.find('ATG', ATG_counter1)
            ATG2 = ref_sequence.find('ATG', ATG_counter2)
            ATG_temp1 = ATG_counter1
            ATG_temp2 = ATG_counter2
            while (ATG1 > 0 or ATG2 > 0) and ATG_limit > 0:
                if (coding_start > 0 and ATG1 > 0):
                    ATG_limit -= 1
                if ATG1 > 0 and ATG2 < 0:
                    ATGs.append([ATG1, -1, ATG1 >= coding_start, seq_previous])
                    ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                elif ATG1 < 0 and ATG2 > 0:
                    ATGs.append([-1, ATG2, ATG2 >= ref_start, seq_previous])
                    ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                elif ATG1-ATG_temp1 == ATG2-ATG_temp2:
                    ATGs.append([ATG1, ATG2, ATG2 >= ref_start, seq_previous])
                    ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                    ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                elif ATG1-ATG_temp1 < ATG2-ATG_temp2:
                    ATGs.append([ATG1, -1, ATG1 >= coding_start, seq_previous])
                    ATG_counter1 = max(ATG_counter1, ATG1 + 1)
                else:
                    ATGs.append([-1, ATG2, ATG2 >= ref_start, seq_previous])
                    ATG_counter2 = max(ATG_counter2, ATG2 + 1)
                ATG1 = sequence.find('ATG', ATG_counter1)
                ATG2 = ref_sequence.find('ATG', ATG_counter2)
            ATG_counter1 = max(0, len(sequence)-2)
            ATG_counter2 = max(0, len(ref_sequence)-2)
            if seq != []:
                seq_previous = seq[2]
            # find transcript-relative coordinates of start codons
            # flatten strings from annotated and reference seqs 
            if seq[1] == 'R':
                if coding_start < 0 and seq[4]*strand + len(seq[0]) > start*strand:
                    coding_start = counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                    ref_start = ref_counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                sequence += seq[0]
                ref_sequence += seq[0]
                counter += len(seq[0])
                ref_counter += len(seq[0])
                continue
            elif seq[2][0][2] == 'D':
                if ref_start < 0 and seq[4]*strand + len(seq[2][0][1]) > start*strand:
                    coding_start = counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                    ref_start = ref_counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                if ((seq[1] == 'G' and include_germline == 1) or 
                    (seq[1] == 'S' and include_somatic == 1)):                  
                    ref_sequence += seq[2][0][1]
                    ref_counter += len(seq[2][0][1])
                continue
            elif seq[2][0][2] == 'I':
                if coding_start < 0 and seq[4]*strand + len(seq[0]) > start*strand:
                    coding_start = counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                    ref_start = ref_counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                sequence += seq[0]
                counter += len(seq[0])
                if ((seq[1] == 'G' and include_germline == 2) or 
                    (seq[1] == 'S' and include_somatic == 2)):                  
                    ref_sequence += seq[0]
                    ref_counter += len(seq[0])
                continue
            elif seq[2][0][2] == 'V':
                if coding_start < 0 and seq[4]*strand + len(seq[0]) > start*strand:
                    coding_start = counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                    ref_start = ref_counter + (start - seq[4] + 1 + 2*self.rev_strand)*strand
                sequence += seq[0]
                counter += len(seq[0])
                if ((seq[1] == 'G' and include_germline == 2) or 
                    (seq[1] == 'S' and include_somatic == 2)):                  
                    ref_sequence += seq[0]
                else:
                    ref_sequence += seq[2][0][1]
                ref_counter += len(seq[0])
                continue
                ##### NEED TO CONSIDER WHAT COULD HAPPEN HERE TO INTRODUCE A NEW
                #### START CODON . . . for instance:
                ## 1 - insert an entirely new codon or piece of one
                ## 2 - Variant could change sequence to create an upstream ATG
                ## 3 - deletion could create new upstream ATG
                ## 4 - modification of original start codon by any of these
                ## THIS NEEDS TO HANDLE NEW CODING START CALC STILL . . .
            # need to process here re: aberrant start codon!!!
            # determine why start codon missing . . . and whether upstream context has changed
            # if no changes upstream and start codon deleted/altered, then 
            # find next downstream to use (even if more are already existing upstream)
            # if novel start introduced upstream, then use that one!
            # changes upstream???
#                    reading_frame = new_start % 3
#                    coding_start += new_start
#                    if reading_frame != 0:
#                        frame_shifts.append([start, -1, 0, -1, []]) 
        print sequence[coding_start:]
        print ref_sequence[ref_start:]
        return []           
        reading_frame = 0
        coordinates = []
        frame_shifts = []
        counter = ref_counter = 0 # hold edited transcript level coordinates
        for seq in annotated_seq:
            # skip sequence fragments that are not to be reported 
            if seq[1] == 'R':
                counter += len(seq[0])
                ref_counter += len(seq[0])
                continue
            elif seq[1] == 'S' and include_somatic < 2:
                if seq[2][0][2] != 'D':
                    counter += len(seq[0])
                    if seq[2][0][2] != 'I':
                        ref_counter += len(seq[0])
                else:
                    ref_counter += len(seq[2][0][1])
                continue
            elif seq[1] == 'G' and include_germline < 2:
                if seq[2][0][2] != 'D':
                    counter += len(seq[0])
                    if seq[2][0][2] != 'I':
                        ref_counter += len(seq[0])
                else:
                    ref_counter += len(seq[2][0][1])
                continue
            # skip sequence fragments that occur prior to start codon 
            # handle cases where variant involves start codon
            if counter < coding_start:
                if seq[2][0][2] == 'D':
                    ref_counter += len(seq[2][0][1])
                    continue
                elif seq[2][0][2] == 'I':
                    if counter + len(seq[0]) > coding_start:
                        coordinates.append([start, seq[4] + len(seq[0])*strand,
                    0, counter + len(seq[0]) - coding_start, seq[2]])
                        if reading_frame != 0:
                            frame_shifts.append([start, -1, 0, -1, []])
                    counter += len(seq[0])
                    continue
                elif seq[2][0][2] == 'V':
                    if counter + len(seq[0]) > coding_start:
                        coordinates.append([start, seq[4] + len(seq[0])*strand,
                    0, counter + len(seq[0]) - coding_start, seq[2]])
                        if reading_frame != 0:
                            frame_shifts.append([start, -1, 0, -1, []])
                    counter += len(seq[0])
                    ref_counter += len(seq[0])
                    continue
                else:
                    # other variant types not handled at this time
                    break                        
            # handle potential frame shifts from indels
            if seq[2][0][2] == 'D' or seq[2][0][2] == 'I':
                if reading_frame == 0:
                    reading_frame = (reading_frame + len(seq[2][0][1])) % 3
                    if reading_frame != 0:
                        frame_shifts.append([seq[2][0][0], -1, counter, -1,seq[2]])
                else:
                    reading_frame = (reading_frame + len(seq[2][0][1])) % 3
                    if reading_frame == 0:
                        # close out all frame_shifts ending in -1
                        for i in range(len(frame_shifts), 0, -1):
                            if frame_shifts[i][1] < 0:
                                frame_shifts[i][1] = seq[4] + len(seq[0])
                                frame_shifts[i][3] = counter + len(seq[0])
                            else:
                                break
                    elif len(seq[2][0][1]) % 3 != 0:
                        frame_shifts.append([seq[2][0][0], -1, counter, -1,seq[2]])
            # log variants                    
            coordinates.append([seq[4], seq[4] + len(seq[0])*strand,
                counter, counter + len(seq[0]), seq[2]])
            if seq[2][0][2] != 'D':
                counter += len(seq[0])
                if seq[2][0][2] != 'I':
                    ref_counter += len(seq[0])
            else:
                ref_counter += len(seq[2][0][1])
        # frame shifts (if they exist) continue to end of transcript
        if reading_frame != 0:
            for i in range(len(frame_shifts), 0, -1):
                if frame_shifts[i][1] < 0:
                    frame_shifts[i][1] = seq[4] + len(seq[0])
                    frame_shifts[i][3] = counter
                else:
                    break

        print sequence
        print sequence[coding_start:]
        print ref_sequence
        print ref_sequence[ref_start:]
        ##### work to be done below here still re: start coordinates and more
        protein = seq_to_peptide(sequence[coding_start:], reverse_strand=False)
        peptide_seqs = kmerize_peptide(seq_to_peptide(sequence[coding_start:], 
            reverse_strand=False), min_size=min_size, max_size=max_size)
        reference_seqs = kmerize_peptide(seq_to_peptide(ref_sequence[ref_start:],
            reverse_strand=False), min_size=min_size, max_size=max_size)
        # get amino acid ranges for kmerization
 ###       for size in range(min_size, max_size + 1):
 ###           epitope_coords = []
 ###           for coords in coordinates:
 ###               # future devel: 
 ###               #    can propogate variant ID here to maintain link to epitope
 ###               epitope_coords.append([max(0, ((coords[0]-start) // 3)-size+1), 
 ###                   min(len(protein), ((coords[1] - start) // 3)+size)])
 ###           for coords in frame_shifts:
 ###               epitope_coords.append([max(0, ((coords[0]-start) // 3)-size+1), 
 ###                   min(len(protein), ((coords[1] - start) // 3)+size)])
 ###           for coords in epitope_coords:
 ###               peptide_seqs += kmerize_peptide(protein[coords[0]:coords[1]], 
 ###                   min_size=size, max_size=size)
        # return list of unique neoepitope sequences
        return list(set(peptide_seqs).difference(reference_seqs))

def gtf_to_cds(gtf_file, dictdir, pickle_it=True):
    """ References cds_dict to get cds bounds for later Bowtie query

        Keys in the dictionary are transcript IDs, while entries are lists of
            relevant CDS/stop codon data
            Data: [chromosome, start, stop, +/- strand]
        Writes cds_dict as a pickled dictionary

        gtf_file: input gtf file to process
        dictdir: path to directory to store pickled dicts

        Return value: dictionary
    """
    cds_dict = collections.defaultdict(list)
    # Parse GTF to obtain CDS/stop codon info
    with open(gtf_file, "r") as f:
        for line in f:
            if line[0] != '#':
                tokens = line.strip().split('\t')
                if tokens[2] in ['exon', 'start_codon', 'stop_codon']: 
                    transcript_id = re.sub(
                                r'.*transcript_id \"([A-Z0-9._]+)\"[;].*', 
                                r'\1', tokens[8]
                                )
                    # Create new dictionary entry for new transcripts
                    cds_dict[transcript_id].append([tokens[0].replace(
                                                                    "chr", ""),
                                                    tokens[2], int(tokens[3]), 
                                                    int(tokens[4]), tokens[6]])
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
            start = cds[2]
            stop = cds[3]
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

def get_transcripts_from_tree(chrom, start, stop, cds_tree):
    """ Uses cds tree to btain transcript IDs from genomic coordinates
            
        chrom: (String) Specify chrom to use for transcript search.
        start: (Int) Specify start position to use for transcript search.
        stop: (Int) Specify ending position to use for transcript search
        cds_tree: (Dict) dictionary of IntervalTree() objects containing
            transcript IDs as function of exon coords indexed by chr/contig ID.
            
        Return value: (set) a set of matching unique transcript IDs.
    """
    transcript_ids = []
    # Interval coordinates are inclusive of start, exclusive of stop
    if chrom not in cds_tree:
        return []
    cds = list(cds_tree[chrom].search(start, stop))
    for cd in cds:
        if cd.data not in transcript_ids:
            transcript_ids.append(cd.data)
    return transcript_ids

if __name__ == '__main__':
    import unittest
    import os
    class TestTranscript(unittest.TestCase):
        """Tests transcript object construction"""
        def setUp(self):
            """Sets up gtf file and creates dictionaries for tests"""
            self.gtf = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'test', 'Chr11.gtf'
                            )
            self.cds = gtf_to_cds(self.gtf, 'NA', pickle_it=False)
            self.ref_prefix = os.path.join(
                            os.path.dirname(
                                    os.path.dirname(
                                            os.path.realpath(__file__)
                                        )
                                ), 'test', 'Chr11.ref'
                        )
            self.reference_index = bowtie_index.BowtieIndexReference(
                                                        self.ref_prefix)
            self.transcript = Transcript(self.reference_index, 
                                        [[str(chrom), 'blah', seq_type,
                                          str(start), str(end), '.', 
                                          strand] for (chrom, seq_type, start, 
                                                        end, strand) in 
                                          self.cds['ENST00000335295.4_1']])
            self.fwd_transcript = Transcript(self.reference_index, 
                                        [[str(chrom), 'blah', seq_type,
                                          str(start), str(end), '.', 
                                          strand] for (chrom, seq_type, start, 
                                                        end, strand) in 
                                          self.cds['ENST00000308020.5_1']])
        def test_transcript_structure(self):
            """Fails if structure of unedited transcript is incorrect"""
            self.assertEqual(len(self.transcript.annotated_seq()), 1)
            self.assertEqual(len(self.transcript.annotated_seq()[0][0]), 628)
            self.assertEqual(self.transcript.annotated_seq()[0][1], 'R')
            self.assertEqual(self.transcript.intervals, [5246692, 5246955, 
                                                         5247805, 5248028,
                                                         5248158, 5248300])
            self.assertEqual(self.transcript.start_codon, 5248248)
            self.assertTrue(self.transcript.rev_strand)
            self.assertEqual(self.transcript.edits, {})
            self.assertEqual(self.transcript.deletion_intervals, [])
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
        def test_irrelevant_edit(self):
            """Fails if edit is made for non-exon position"""
            self.transcript.edit('G', 5248155)
            relevant_edits = self.transcript.expressed_edits()
            self.assertEqual(self.transcript.edits[5248154], [('G', 'V', 'S', 
                                                              (5248155, 'G', 
                                                                'V'), None)])
            self.assertEqual(relevant_edits[0], {})
            self.assertEqual(relevant_edits[1], [(5246692, 'R', '', None),
                                                 (5246955, 'R', '', None),
                                                 (5247805, 'R', '', None), 
                                                 (5248028, 'R', '', None), 
                                                 (5248158, 'R', '', None), 
                                                 (5248300, 'R', '', None)])
            self.assertEqual(len(self.transcript.annotated_seq()), 1)
        def test_relevant_edit(self):
            """Fails if edit is not made for position within exon"""
            self.transcript.edit('A', 5248299)
            relevant_edits = self.transcript.expressed_edits()
            self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 'S',
                                                              (5248299, 'A', 
                                                                'V'), None)])
            self.assertEqual(relevant_edits[0][5248298], [('A', 'V', 'S',
                                                          (5248299, 'A', 
                                                                'V'), None)])
            self.assertEqual(relevant_edits[1], [(5246692, 'R', '', None),
                                                 (5246955, 'R', '', None),
                                                 (5247805, 'R', '', None), 
                                                 (5248028, 'R', '', None), 
                                                 (5248158, 'R', '', None), 
                                                 (5248300, 'R', '', None)])
            self.assertEqual(len(self.transcript.annotated_seq()), 3)
            self.assertEqual(self.transcript.annotated_seq()[1][0], 'T')
            self.assertEqual(self.transcript.annotated_seq()[1][1], 'S')
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
                                                                   (5248299, 
                                                                    'A', 'V'),
                                                                    None)])
            self.assertEqual(self.transcript.last_deletion_intervals,
                                [(5246692, 5246695, 'S', (5246694, 'TTG', 
                                                           'D'), None)])
        def test_reset_to_save_point(self):
            """Fails if new edit not erased or old edits not retained"""
            self.transcript.edit('A', 5248299)
            self.transcript.edit(3, 5246694, mutation_type='D')
            self.transcript.save()
            self.transcript.edit('G', 5248165)
            self.assertEqual(self.transcript.edits[5248164], [('G', 'V', 'S',
                                                              (5248165, 'G',
                                                                'V'), None)])
            self.transcript.reset(reference=False)
            self.assertNotIn(2182387, self.transcript.edits)
            self.assertEqual(self.transcript.last_edits[5248298], [('A', 'V', 
                                                                    'S', 
                                                                   (5248299, 
                                                                    'A', 'V'),
                                                                    None)])
            self.assertEqual(self.transcript.last_deletion_intervals,
                                [(5246692, 5246695, 'S', (5246694, 'TTG', 
                                                           'D'), None)])
            self.assertNotEqual(self.transcript.edits, {})
        def test_SNV_seq(self):
            """Fails if SNV is edited incorrectly"""
            self.transcript.edit('A', 5248299)
            self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 
                                                                    'S', 
                                                                   (5248299, 
                                                                    'A', 'V'),
                                                                    None)])
            self.assertEqual(self.transcript.deletion_intervals, [])
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(seq[0], ('AC', 'R', [''], [], 5248301))
            self.assertEqual(seq[1], ('T', 'S', [(5248299, 'A', 'V')], [], 
                                        5248299))
            self.assertEqual(len(seq[2][0]), 625)
        def test_inside_insertion(self):
            """Fails if indel within exon is inserted incorrectly"""
            self.transcript.edit('Q', 5248165, mutation_type='I')
            self.assertEqual(self.transcript.edits[5248164], [('Q', 'I', 'S',
                                                                (5248165, 'Q',
                                                                'I'), None)])
            self.assertEqual(self.transcript.deletion_intervals, [])
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(seq[1], ('Q', 'S', [(5248165, 'Q', 'I')], [], 
                                        5248165))
            self.assertEqual(len(seq[0][0]), 136)
            self.assertEqual(len(seq[2][0]), 492)
        def test_adjacent_indel(self):
            """Fails if indel right before exon is inserted incorrectly"""
            self.transcript.edit('Q', 5248029, mutation_type='I')
            self.assertEqual(self.transcript.edits[5248028], [('Q', 'I', 'S',
                                                                (5248029, 'Q',
                                                                'I'), None)])
            self.assertEqual(self.transcript.deletion_intervals, [])
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(seq[1], ('Q', 'S', [(5248029, 'Q', 'I')], [], 
                                        5248029)) 
            self.assertEqual(len(seq[0][0]), 142)
            self.assertEqual(len(seq[2][0]), 486)
        def test_inside_deletion(self):
            """Fails if deletion completely within an exon is made improperly"""
            self.transcript.edit(3, 5246700, mutation_type='D')
            self.assertEqual(self.transcript.edits, {})
            self.assertEqual(self.transcript.deletion_intervals, [(5246698, 
                                                                   5246701, 'S',
                                                                   (5246700, 
                                                                    'TGA', 'D'), 
                                                                   None)])
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(seq[1], ('', 'S', [(5246700, 'TGA', 'D')], [], 
                                        5246700))
            self.assertEqual(seq[2], ('TTGCAA', 'R', [''], [], 5246700))
            self.assertEqual(len(seq[0][0]), 619)
        def test_overlapping_deletion(self):
            """Fails if deletion overlapping a junction is incorrect"""
            self.transcript.edit(10, 5246950, mutation_type='D')
            self.assertEqual(self.transcript.edits, {})
            self.assertEqual(self.transcript.deletion_intervals, [(5246948,
                                                                   5246958, 'S',
                                                                   (5246950,
                                                                    'CCAGGAGCTG',
                                                                    'D'),
                                                                   None)])
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(seq[1], ('', 'S', [(5246950,'CCAGGAGCTG', 'D')], 
                                        [], 5246950))
            self.assertEqual(len(seq[0][0]), 365)
            self.assertEqual(len(seq[2][0]), 256)
        def test_spanning_deletion(self):
            self.transcript.edit(137, 5248025, mutation_type="D")
            self.assertEqual(self.transcript.edits, {})
            self.assertEqual(self.transcript.deletion_intervals[0][0:3],
                                                    (5248023, 5248160, 'S'))
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 3)
            self.assertEqual(len(seq[1][2][0][1]), 137)
            self.assertEqual(seq[1][0:2], ('', 'S'))
            self.assertEqual(len(seq[0][0]), 140)
            self.assertEqual(len(seq[2][0]), 481)
        def test_compound_variants(self):
            """Fails if transcript with multiple variant types is incorrect"""
            self.transcript.edit(137, 5248025, mutation_type='D')
            self.transcript.edit('Q', 5248165, mutation_type='I')
            self.transcript.edit('A', 5248299)
            self.assertEqual(len(self.transcript.edits.keys()), 2)
            self.assertEqual(self.transcript.edits[5248298], [('A', 'V', 
                                                                    'S', 
                                                                   (5248299, 
                                                                    'A', 'V'),
                                                                    None)]) 
            self.assertEqual(self.transcript.edits[5248164], [('Q', 'I', 'S',
                                                                (5248165, 'Q',
                                                                'I'), None)])
            self.assertEqual(self.transcript.deletion_intervals[0][0:3],
                                                    (5248023, 5248160, 'S'))
            seq = self.transcript.annotated_seq()
            self.assertEqual(len(seq), 7)
            self.assertEqual(seq[0], ('AC', 'R', [''], [], 5248301))
            self.assertEqual(seq[1], ('T', 'S', [(5248299, 'A', 'V')], [], 
                                        5248299))
            self.assertEqual(len(seq[2][0]), 133)
            self.assertEqual(seq[3], ('Q', 'S', [(5248165, 'Q', 'I')], [], 
                                        5248165))
            self.assertEqual(seq[4], ('GGGC', 'R', [''], [], 5248165))
            self.assertEqual(len(seq[6][0]), 481)
    unittest.main()
