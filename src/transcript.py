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

    def edit(self, seq, pos, mutation_type='V', mutation_class='S'):
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
                            seq)
                    )
            else:
                self.deletion_intervals.append(
                        (pos - 2, pos + deletion_size - 2, mutation_class, 
                            self.bowtie_reference_index.get_stretch(
                                    self.chrom, pos - 1, 
                                    pos + deletion_size + 2 - pos - 1
                                ))
                    )
        elif mutation_type == 'V' or mutation_type == 'I':
            self.edits[pos - 1].append((seq, mutation_type, mutation_class))
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
                                   sorted_deletion_intervals[0][2]),
                                  (sorted_deletion_intervals[0][1],
                                   sorted_deletion_intervals[0][2])]
            for i in xrange(1, len(sorted_deletion_intervals)):
                if (sorted_deletion_intervals[i][0]
                    <= deletion_intervals[-1][0]):
                    deletion_intervals[-2] = min(deletion_intervals[-2],
                                            (sorted_deletion_intervals[i][0],
                                             sorted_deletion_intervals[i][2]),
                                            key=itemgetter(0))
                    deletion_intervals[-1] = max(deletion_intervals[-1],
                                            (sorted_deletion_intervals[i][1],
                                             sorted_deletion_intervals[i][2]),
                                            key=itemgetter(0))
                else:
                    deletion_intervals.extend(
                            [(sorted_deletion_intervals[i][0],
                                sorted_deletion_intervals[i][2]),
                             (sorted_deletion_intervals[i][1],
                                sorted_deletion_intervals[i][2])]
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
                else:
                    assert end_index > start_index
                    if start_index % 2:
                        pos = deletion_intervals[i]
                    else:
                        pos = (intervals[start_index], 'R')
                        start_index += 1
                    # deletion_intervals[i] becomes a new end
                    relevant_deletion_intervals.extend(
                            [pos, (intervals[start_index], 'R')]
                        )
                    if end_index % 2:
                        end_pos = deletion_intervals[i+1]
                        relevant_deletion_intervals.extend(
                            [(intervals[i], 'R') for i in
                             xrange(start_index + 1, end_index)]
                        )
                    else:
                        end_pos = (intervals[end_index - 1], 'R')
                        relevant_deletion_intervals.extend(
                                [(intervals[i], 'R') for i in
                                 xrange(start_index, end_index)]
                            )
                    relevant_deletion_intervals.append(end_pos)
        intervals = sorted([(interval, 'R') for interval in intervals]
                            + relevant_deletion_intervals)
        edits = collections.defaultdict(list)
        for pos in self.edits:
            # Add edit if and only if it's in one of the CDSes
            start_index = custom_bisect_left(intervals, pos)
            for edit in self.edits[pos]:
                if (include_somatic and edit[2] == 'S'
                        or include_germline and edit[2] == 'G'):
                    if edit[1] == 'V':
                        if start_index % 2:
                            # Add edit if and only if it lies within bounds
                            edits[pos].append(edit)
                    elif edit[1] == 'I':
                        if start_index % 2 or pos == intervals[start_index][0]:
                            # An insertion is valid before or after a block
                            edits[pos].append(edit)
        return (edits, intervals)

    def save(self):
        """ Creates save point for edits.

            No return value.
        """
        self.last_edits = copy.copy(self.edits)
        self.last_deletion_intervals = copy.copy(self.deletion_intervals)

    def _seq_append(self, seq_list, seq, mutation_class):
        """ Appends mutation to seq_list, merging successive mutations.

            seq_list: list of tuples (sequence, type) where type is one
                of R, G, or S (for respectively reference, germline edit, or
                somatic edit). Empty sequence means there was a deletion.
            seq: seq to add
            mutation_class: S for somatic, G for germline, R for reference

            No return value; seq_list is merely updated.
        """
        try:
            condition = seq_list[-1][-1] == mutation_class
        except IndexError:
            # Add first item in seq_list
            assert not seq_list
            if seq or mutation_class != 'R':
                seq_list.append((seq, mutation_class))
            return
        if condition:
            seq_list[-1] = (seq_list[-1][0] + seq, mutation_class)
        elif seq or mutation_class != 'R':
            seq_list.append((seq, mutation_class))

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

            Return value: list of triples (sequence, variant, type) where 
                variant is one of V, I, or D (for SNV, insertion, or deletion, 
                respectively) and type is one of R, G, or S (for reference, 
                germline edit, or somatic edit, respectively). Sequence will be 
                deletion size (in bp) for deletion variants.
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
                if intervals[i][0] in edits:
                    assert (len(edits[intervals[i][0]]) == 1
                                and edits[intervals[i][0]][0][1] == 'I')
                    if i:
                        new_edits[
                            intervals[i-1][0]] = new_edits[intervals[i][0]]
                        del new_edits[intervals[i][0]]
                    else:
                        intervals = [(-1, 'R'), (-1, 'R')] + intervals
                        # Have to add 2 because we modified intervals above
                        i += 2
                        new_edits[-1] = new_edits[intervals[i][0]]
                        del new_edits[intervals[i][0]]
                i += 2
            seqs = []
            for i in xrange(0, len(intervals), 2):
                seqs.append(
                        self.bowtie_reference_index.get_stretch(
                                self.chrom, intervals[i][0] + 1,
                                intervals[i + 1][0] -
                                intervals[i][0]
                            )
                    )
            # Now build sequence in order of increasing edit position
            i = 1
            pos_group, final_seq = [], []
            for pos in (sorted(new_edits.keys()) + [(self.intervals[-1] + 1,
                                                        'R')]):
                if pos > intervals[i][0]:
                    last_index, last_pos = 0, intervals[i-1][0] + 1
                    for pos_to_add in pos_group:
                        fill = pos_to_add - last_pos
                        if intervals[i-1][1] != 'R':
                            self._seq_append(final_seq, '', intervals[i-1][1])
                        self._seq_append(final_seq, seqs[(i-1)/2][
                                            last_index:last_index + fill
                                        ], 'R')
                        if intervals[i][1] != 'R':
                            self._seq_append(final_seq, '', intervals[i][1])
                        # If no edits, snv is reference and no insertion
                        try:
                            snv = (seqs[(i-1)/2][last_index + fill], 'R')
                        except IndexError:
                            '''Should happen only for insertions at beginning
                            of sequence.'''
                            assert (i - 1) / 2 == 0 and not seqs[0]
                            snv = ('', 'R')
                        insertion = ('', 'R')
                        for edit in new_edits[pos_to_add]:
                            if edit[1] == 'V':
                                snv = (edit[0], edit[2])
                            else:
                                assert edit[1] == 'I'
                                insertion = (edit[0], edit[2])
                        self._seq_append(final_seq, *snv)
                        self._seq_append(final_seq, *insertion)
                        last_index += fill + 1
                        last_pos += fill + 1
                    if intervals[i-1][1] != 'R':
                        self._seq_append(final_seq, '', intervals[i-1][1])
                    self._seq_append(
                            final_seq, seqs[(i-1)/2][last_index:], 'R'
                        )
                    if intervals[i][1] != 'R':
                        self._seq_append(final_seq, '', intervals[i][1])
                    i += 2
                    try:
                        while pos > intervals[i]:
                            if intervals[i-1][1] != 'R':
                                self._seq_append(
                                        final_seq, '', intervals[i-1][1]
                                    )
                            self._seq_append(final_seq, seqs[(i-1)/2], 'R')
                            if intervals[i][1] != 'R':
                                self._seq_append(
                                        final_seq, '', intervals[i][1]
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
                return [(seq[::-1].translate(revcomp_translation_table),
                            mutation_class)
                            for seq, mutation_class in final_seq][::-1]
            return final_seq
        raise NotImplementedError(
            'Retrieving sequence with transcript coordinates not '
            'yet fully supported.'
        )

    def peptides(self, size=9, somatic=2, germline=1):
        """ Retrieves list of predicted peptide fragments from transcript that 
            include one or more variants.

            size: peptide length (specified as # of amino acids)
            somatic: 0 to omit consideration of somatic variants, 1 to identify 
            somatic variants but not explicitly report peptides containing every
            somatic variant, 2 to identify and report all peptides containing 
            consequences of somatic variants
            germline: 0 to omit consideration of germline variants, 1 to 
            identify germline variants but not explicitly report peptides 
            containing every germline variant, 2 to identify and report all 
            peptides containing consequences of germline variants

            Return value: list of peptides of desired length.
        """
        if size < 2: return []
        annotated_seq = self.annotated_seq(include_somatic=somatic != 0, 
            include_germline=germline != 0)
        coordinates = []
        counter = 0 # hold transcript level coordinates
        frame_shifts = []
        sequence = '' # hold flattened nucleotide sequence
        # extract nucleotide sequence from annotated_seq
        for seq in annotated_seq:
            record = (seq[2] == 'S' and somatic >= 2) or (seq[2] == 'G' 
                and germline >= 2)
            if seq[1] == 'D' and record:
                coordinates.append((counter, 0))
        # locate position of start codon (first ATG in sequence)
        start = sequence.find("ATG")
        if start < 0: return []
        reading_frame = (start - self.start_codon) % 3
        if reading_frame != 0:
            frame_shifts.append((start, start))
        for seq in annotated_seq:
            # skip sequence fragments that occur prior to start codon 
            if seq[1] != 'D' and counter + len(seq[0]) <= start:
                counter += len(seq[0])
                continue
            elif seq[1] == 'D' and counter < start:
                continue
            # skip sequence fragments that are not to be reported 
            if (seq[2] == 'R' or (seq[2] == 'S' and somatic < 2) or 
                (seq[2] == 'G' and germline >= 2)):
                if seq[1] != 'D':
                    counter += len(seq[0])
                continue
            # handle unique case where variant precedes but includes start codon
            if counter < start:
                coordinates.append(start, counter + len(seq[0]))
                if seq[1] == 'I' and reading_frame == 0:
                    reading_frame = (reading_frame + len(seq[0])) % 3
                    if reading_frame != 0:
                        frame_shifts.append((counter, counter))
                elif seq[1] == 'I':
                    reading_frame = (reading_frame + len(seq[0])) % 3
                    if reading_frame == 0:
                        frame_shifts[-1][1] = counter + len(seq[0]) 
                counter += len(seq[0])                  #
                continue
            # handle potential frame shifts from indels
            if seq[1] == 'D' or seq[1] == 'I':
                if reading_frame == 0:
                    reading_frame = (reading_frame + seq[0]) % 3
                    if reading_frame != 0:
                        frame_shifts.append((counter, counter))
                else:
                    reading_frame = (reading_frame + seq[0]) % 3
                    if reading_frame == 0:
                        frame_shifts[-1][1] = counter
            # log variants                    
            if seq[1] == 'D':
                coordinates.append((counter, 0))
            else:
                coordinates.append((counter, len(seq[0])))
                counter += len(seq[0])
            continue
        # frame shift (if it exists) continues to end of transcript
        if reading_frame != 0:
            frame_shifts[-1][1] = counter
        protein = seq_to_peptide(sequence[start:])
        # for each variant and any areas of different reading frame, do the windows around there

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
    cds_dict = {}
    # Parse GTF to obtain CDS/stop codon info
    with open(gtf_file, "r") as f:
        for line in f:
            if line[0] != '#':
                tokens = line.strip().split('\t')
                if tokens[2] == "exon":
                    transcript_id = re.sub(
                                r'.*transcript_id \"([A-Z0-9._]+)\"[;].*', 
                                r'\1', tokens[8]
                                )
                    # Create new dictionary entry for new transcripts
                    if transcript_id not in cds_dict:
                        cds_dict[transcript_id] = [[tokens[0].replace(
                                                                    "chr", ""), 
                                                        int(tokens[3]), 
                                                        int(tokens[4]), 
                                                        tokens[6]]]
                    else:
                        cds_dict[transcript_id].append([tokens[0].replace(
                                                                    "chr", ""), 
                                                            int(tokens[3]), 
                                                            int(tokens[4]), 
                                                            tokens[6]])
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

def get_transcripts_from_tree(chrom, start, stop, cds_tree):
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
    cds = cds_tree[chrom].search(start, stop)
    for cd in cds:
        if cd.data not in transcript_ids:
            transcript_ids.add(cd.data)
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
                                    ), 'test', 'Ychrom.gtf'
                            )
            self.cds = gtf_to_cds(self.gtf, "NA", pickle_it=False)
            self.ref_prefix = os.path.join(
                            os.path.dirname(
                                    os.path.dirname(
                                            os.path.realpath(__file__)
                                        )
                                ), 'test', 'Ychrom.ref'
                        )
            self.reference_index = bowtie_index.BowtieIndexReference(
                                                        self.ref_prefix)
            self.CDS_lines = self.cds['ENST00000421783.1_2']
            self.transcript = Transcript(self.reference_index, 
                                        [[str(chrom), 'blah', 'blah',
                                          str(start), str(end), '.', 
                                          strand] for (chrom, start, 
                                                        end, strand) in 
                                          self.CDS_lines])
        def test_transcript_structure(self):
            """Fails if structure of unedited transcript is incorrect"""
            self.assertEqual(len(self.transcript.annotated_seq()), 1)
            self.assertEqual(len(self.transcript.annotated_seq()[0][0]), 464)
            self.assertEqual(self.transcript.annotated_seq()[0][1], 'R')
            self.assertEqual(self.transcript.intervals, [2181011, 2181101, 
                                                         2182013, 2182387])
            self.assertFalse(self.transcript.rev_strand)
            self.assertEqual(self.transcript.edits, {})
            self.assertEqual(self.transcript.deletion_intervals, [])
        def test_irrelevant_edit(self):
            """Fails if edit is made for non-exon position"""
            self.transcript.edit('G', 2181009)
            relevant_edits = self.transcript.expressed_edits()
            self.assertEqual(self.transcript.edits[2181008], [('G', 'V', 'S')])
            self.assertEqual(relevant_edits[0], {})
            self.assertEqual(relevant_edits[1], [(2181011, 'R'), 
                                                 (2181101, 'R'),
                                                 (2182013, 'R'),
                                                 (2182387, 'R')])


    '''
        def test_relevant_edit(self):
            """Fails if edit is not made for CDS position"""
            self.transcript.edit("G", 34)
            relevant_edits = self.transcript.expressed_edits()
            self.assertEqual(self.transcript.edits[33], [("G", "V", "S")])
            self.assertEqual(relevant_edits[0][33], [("G", "V", "S")])
            self.assertEqual(relevant_edits[1], [(29, "R"), (74, "R"), 
                                                    (74, "R"), (77, "R")])
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
            self.assertEqual(self.transcript.last_edits[33], [("G", "V", 
                                                                "S")])
            self.assertEqual(self.transcript.last_deletion_intervals,
                                [(58, 61, "S")])
        def test_reset_to_save_point(self):
            """Fails if new edit not erased or old edits not retained"""
            self.transcript.edit("G", 34)
            self.transcript.edit("CCC", 60, mutation_type="D")
            self.transcript.save()
            self.transcript.edit("C", 36)
            self.assertEqual(self.transcript.edits[35], [("C", "V", "S")])
            self.transcript.reset(reference=False)
            self.assertNotIn(35, self.transcript.edits)
            self.assertEqual(self.transcript.last_edits[33], [("G", "V", 
                                                                "S")])
            self.assertEqual(self.transcript.last_deletion_intervals,
                                [(58, 61, "S")])
            self.assertNotEqual(self.transcript.edits, {})
        def test_SNV_seq(self):
            """Fails if SNV is edited incorrectly"""
            self.transcript.edit("G", 31)
            seq1 = self.transcript.annotated_seq()
            seq2 = self.transcript.annotated_seq(31, 36)
            self.assertEqual(seq1, [('G', 'S'), 
                ('TGCCCGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG', 'R')])
            self.assertEqual(seq2, [('G', 'S'), ('TGCCC', 'R')])
        def test_inside_indel(self):
            """Fails if indel within CDS is inserted incorrectly"""
            self.transcript.edit("Q", 35, mutation_type="I")
            self.assertEqual(self.transcript.edits[34], [("Q", "I", "S")])
            seq1 = self.transcript.annotated_seq()
            seq2 = self.transcript.annotated_seq(31, 36)
            self.assertEqual(seq1, [('ATGCC', 'R'), ('Q', 'S'), 
                    ('CGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG', 'R')])
            self.assertEqual(seq2, [('ATGCC', 'R'), ('Q', 'S'), ('C', 'R')])
        def test_adjacent_indel(self):
            """Fails if indel right before CDS is inserted incorrectly"""
            self.transcript.edit("Q", 30, mutation_type="I")
            self.assertEqual(self.transcript.edits[29], [("Q", "I", "S")])
            seq1 = self.transcript.annotated_seq()
            seq2 = self.transcript.annotated_seq(30, 36)
            ## DEBUG THESE
            self.assertEqual(seq1, [('Q', 'S'), 
                ('ATGCCCGTGCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG', 'R')])
            self.assertEqual(seq2, [('Q', 'S'), ('ATGCCC', 'R')])
        def compound_variants(self):
            self.transcript.edit("Q", 30, mutation_type ="I")
            self.transcript.edit("T", 33, mutation_type ="V")
            self.transcript.edit("JJJ", 40, mutation_type ="I")
            self.transcript.edit(1, 45, mutation_type ="D")
            seq1 = self.transcript.annotated_seq()
            seq2 = self.transcript.annotated_seq(30, 36)
            self.assertEqual(seq1, [('Q', 'S'), ('AT', 'R'), ('T', 'S'), 
                                    ('CCCGTGC', 'R'),  ('JJJ', 'S'), 
                                    ('CGAA', 'R'), ('', 'S'), 
                                    ('TCGTGTCCCCGCTACAATGCCCGTGCCGATTTG', 
                                        'R')])
            self.assertEqual(seq2, [('Q', 'S'), ('AT', 'R'), ('T', 'S'), 
                                                            ('CC', 'R')])
        def test_deletion(self):
            """Fails if deletion is made incorrectly"""
            self.transcript.edit(5, 34, mutation_type="D")
            self.assertEqual(self.transcript.deletion_intervals, [(32, 37, 
                                                                    "S")])
            seq1 = self.transcript.annotated_seq()
            seq2 = self.transcript.annotated_seq(31, 36)
            self.assertEqual(seq1, [('ATG', 'R'), ('', 'S'), 
                        ('GCCGAATTCGTGTCCCCGCTACAATGCCCGTGCCGATTTG', 'R')])
            ### DEBUG THIS ONE
            self.assertEqual(seq2, [('ATG', 'R'), ('', 'S'), ('GCC', 'R')])
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
    '''
    unittest.main()
