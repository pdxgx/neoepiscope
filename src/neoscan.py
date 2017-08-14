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


def kmer(mutation_posits, size_min, size_max, normal_aa, mutated_aa = ""):
    ''' Generates kmers of size size_min to size_max (e.g. 8-11) from amino acid sequences.
        mutation_posits: (List) Positions of mutations on chromosome
        normal_aa: (String) Amino acid sequence of Wild Type
        mutated_aa: (String) Amino acid sequence of Mutant Type
        Return value: (List) List of tuples pairing mismatching kmers for a given epitope.
    '''
    if (len(mutated_aa) == 0):
        mutated_aa = normal_aa
    kmer_list = list()
    #Loop through window sizes
    for ksize in range(size_min, size_max+1):
        for startIndex in range(len(mutated_aa)-ksize+1):
            kmer_list.append((normal_aa[startIndex:startIndex+ksize], mutated_aa[startIndex:startIndex+ksize]))
    final_list = list()
    for WT,MT in kmer_list:
        if (WT != MT and 'X' not in MT):
            final_list.append((WT, MT))
    return final_list

def get_cds(transcript_id, mutation_pos_list, seq_length_left, 
              seq_length_right, ordered_cds_dict, mutation_dict):
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
    bounds_set = set()
    if transcript_id not in ordered_cds_dict:
        print('failed on 1')
        print(transcript_id)
        print(ordered_cds_dict)
        return [], mutation_dict, bounds_set
    pos_in_codon = 2 - (seq_length_right%3)
    cds_list = ordered_cds_dict[transcript_id]
    mutation_pos = -1
    #Don't want to check rightmost since seq. queries based off of it.
    if len(mutation_pos_list) >= 2:
        removal_list = []
        shift = mutation_pos_list[0][0] - min(mutation_dict)
        key_list = list(mutation_dict.keys())
        #Remove all mutations outside of cds boundaries.
        for index in range(len(mutation_pos_list)):
            lower_cds_index = 2*bisect.bisect(cds_list[::2], mutation_pos_list[index][0])-2
            upper_cds_index = lower_cds_index+1
            if(lower_cds_index < 0 or 
               cds_list[upper_cds_index] < mutation_pos_list[index][0]):
                #Delete at the current index
                try:
                    if mutation_dict[key_list[index]] == '':
                        try:
                            count = 1
                            while mutation_dict[key_list[index+count]] == '':
                                count += 1
                                removal_list.append(index+count)
                                del mutation_dict[key_list[index+count]]
                        except IndexError:
                            pass
                    del mutation_dict[key_list[index]]
                    removal_list.append(index)
                except IndexError:
                    continue
        for index in range(len(removal_list)-1, -1, -1):
            #print("made edits to mutation pos list")
            mutation_pos_list.pop(removal_list[index])
    #Loop again, this time from right & correcting seq queries.
    for index in range(len(mutation_pos_list)-1, -1, -1):
        mutation = mutation_pos_list[index][0]
        middle_cds_index = 2*bisect.bisect(cds_list[::2], mutation)-2
        #If the middle_cds_index is past the last boundary, move it to the last
        if middle_cds_index > len(cds_list)-1:
            middle_cds_index -= 2
        #If the biggest pogsition is smaller than the smallest bound, return []
        if middle_cds_index < 0:
            print('failed on 2')
            return [], mutation_dict, bounds_set
        curr_left_index = middle_cds_index
        curr_right_index = middle_cds_index+1 #cds boundary end indexes
        #Increase by one to ensure mutation_pos_list is collected into boundary
        curr_pos_left = mutation + 1
        curr_pos_right = mutation #Actual number in chromosome
        #If the mutation is not on in cds bounds, return [].
        if(mutation <= cds_list[curr_right_index] and 
           mutation >= cds_list[curr_left_index]):
            mutation_pos = mutation
            if index != len(mutation_pos_list)-1:
                #shift is the current mutation's position in the codon.
                new_pos_in_codon = (mutation_pos_list[-1][0]
                                    - pos_in_codon-mutation_pos_list[index][0]) % 3
                seq_length_right = (_size_max-1)*3 + new_pos_in_codon
                seq_length_left -= (mutation_pos_list[-1][0] 
                                    - mutation_pos_list[index][0])
            break
    if(mutation_pos == -1):
        print('failed on 3')
        return [], mutation_dict, bounds_set
    #Increase the seq length by 1 to account for mutation_pos_list collection
    seq_length_left += 1
    total_seq_length = seq_length_right + seq_length_left
    original_length_left = seq_length_left
    nucleotide_index_list = []
    #Loop left until receive all queried left-side bases and mutation base.
    while(len(nucleotide_index_list) == 0 or 
          sum([index[1] for index in nucleotide_index_list]) 
          < (original_length_left)):
        if curr_pos_left-cds_list[curr_left_index]+1 >= seq_length_left:
            if curr_pos_left != mutation_pos+1:
                nucleotide_index_list.append((curr_pos_left-seq_length_left+1,
                                          seq_length_left))
            else:

                nucleotide_index_list.append((curr_pos_left-seq_length_left,
                                          seq_length_left))
            bounds_set.add((cds_list[curr_left_index], cds_list[curr_left_index+1]))
            seq_length_left = 0
        else:
            if curr_pos_left != mutation_pos+1:
                curr_pos_left += 1
            nucleotide_index_list.append((cds_list[curr_left_index],
                                    curr_pos_left-cds_list[curr_left_index]))
            bounds_set.add((cds_list[curr_left_index], cds_list[curr_left_index+1]))
            seq_length_left -= (curr_pos_left-cds_list[curr_left_index])
            curr_pos_left = cds_list[curr_left_index-1]
            curr_left_index -= 2
            if curr_left_index < 0:
                #print("Exceeded all possible cds boundaries!")
                #Changed total_seq_length for comparison in next while loop.
                total_seq_length = (original_length_left
                                      - seq_length_left
                                      + seq_length_right)
                break
    #Reverse list to get tuples in order
    nucleotide_index_list = list(reversed(nucleotide_index_list))
    while(len(nucleotide_index_list) == 0 or 
              sum([index[1] for index in nucleotide_index_list]) 
              < (total_seq_length)):
        if cds_list[curr_right_index] >= curr_pos_right + seq_length_right:
            if curr_pos_right == mutation_pos:
                nucleotide_index_list.append((curr_pos_right+1,
                                              seq_length_right))
            else:
                nucleotide_index_list.append((curr_pos_right,
                                              seq_length_right))
            bounds_set.add((cds_list[curr_right_index-1], cds_list[curr_right_index]))
            seq_length_right = 0
        else:
            try:
                if curr_pos_right == mutation:
                    curr_pos_right += 1
                nucleotide_index_list.append((curr_pos_right,
                                              cds_list[curr_right_index]
                                              - curr_pos_right + 1))
                bounds_set.add((cds_list[curr_right_index-1], cds_list[curr_right_index]))
                seq_length_right -= (cds_list[curr_right_index]-curr_pos_right+1)
                curr_pos_right = cds_list[curr_right_index+1]
                curr_right_index += 2
            except IndexError:
                #print("Exceeded all possible cds boundaries!")
                break
    return nucleotide_index_list, mutation_dict, bounds_set

def get_introns(start_index, stop_index, cds_list):
    total_introns = 0
    index = start_index
    while index < stop_index:
        total_introns += cds_list[index+1]-cds_list[index]-1
        index += 2
    return total_introns

def find_stop(query_st, trans_id, line_count, cds_dict, chrom, reference_index, mutation_locs, reverse):
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
    until_stop = ""
    start = query_st
    stop_found = False
    cds_end = False
    (l_query, r_query) = (0, 3*_size_max)
    if reverse:
        (l_query, r_query) = (r_query, l_query)
    while(stop_found == False):
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(start,line_count)], l_query, r_query, cds_dict, mutation_locs)
        extra_cods = ""
        if reverse:
            print('exons', exon_list)
            if len(exon_list) != 0:
                exon_list = exon_list[:-1] + [(exon_list[-1][0], exon_list[-1][1]-1)]
        else:
            if len(exon_list) != 0:
                if len(exon_list) != 1:
                    exon_list = [(exon_list[0][0], exon_list[0][1]-1)] + exon_list[1:]
                else:
                    exon_list = [(exon_list[0][0], exon_list[0][1]-1)]
        for bound_start, bound_stretch in exon_list:
            extra_cods += get_seq(chrom, bound_start, bound_stretch, reference_index)
        if reverse:
        #    extra_cods = extra_cods[:34]
            last_start = start
            start = exon_list[0][0]
        #    print "extra_cods", start, query_st, l_query, extra_cods
        else:
            start = exon_list[-1][0] + exon_list[-1][1] - 1
        print 'extra cods', extra_cods
        if len(extra_cods) != 3*_size_max:
            print 'triggered', len(extra_cods)
            cds_end = True
        count = 0
        while(count<len(extra_cods)):
            if reverse:
                new_codon = extra_cods[len(extra_cods)-3-count:len(extra_cods)-count]
                count += 3
                amino_acid = seq_to_peptide(new_codon, reverse_strand=True)
            else:
                new_codon = extra_cods[count: count+3]
                count += 3
                amino_acid = seq_to_peptide(new_codon)
            if(amino_acid==""):
                print 317, "found stop", new_codon, extra_cods, extra_cods[len(extra_cods)-3-count:len(extra_cods)-count], count
                stop_found = True
                break
            if reverse:
                until_stop = new_codon + until_stop
            else:
                until_stop += new_codon
        if cds_end:
            break
    if reverse:
        #It might be the case that we need to add whatever we had pulled from 3-get_cds%3 to our splice. Why tho
        #if cds_end:
        #    (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(last_start,line_count)], count, 0, cds_dict, mutation_locs)
        #else:
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(last_start,line_count)], count, 0, cds_dict, mutation_locs)
        print "sequence ending position is this", exon_list[0][0], extra_cods, last_start, count, get_seq(chrom, exon_list[0][0], 9, reference_index)
        return until_stop, exon_list[0][0] #@TODO check whether the +3 is correct
    else:
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(start,line_count)], 0, count, cds_dict, mutation_locs)
        return until_stop, exon_list[-1][0]+exon_list[-1][1]-1
    #return until_stop

def get_seq(chrom, start, splice_length, reference_index):
    chr_name = "chr" + chrom #proper
    start -= 1 #adjust for 0-based bowtie queries
    try:
        seq = reference_index.get_stretch(chr_name, start, splice_length)
    except KeyError:
        return False
    return seq

def make_mutation_seq(orig_seq, mutation_locs, reverse):
    mutation_seq = ""
    if(reverse):
        orig_seq = orig_seq[::1]
        for ind in range(len(orig_seq)):
            if ind in mutation_locs:
                mutation_seq = mutation_locs[ind]+mutation_seq
            else:
                mutation_seq = orig_seq[ind] + mutation_seq 
        return mutation_seq[::-1]
    else:
        for ind in range(len(orig_seq)):
            if ind in mutation_locs:
                mutation_seq += mutation_locs[ind]
            else:
                mutation_seq += orig_seq[ind]
        return mutation_seq

def find_seq_and_kmer(cds_list, last_chrom, reference_index, mutation_locs,
                      orf_dict, trans_id, mutation_posits, direct):
    #if bam exists:
    if args.bam:
        seq_start = cds_list[0][0]
        (last_start, last_length) = cds_list.pop()
        seq_end = last_start+last_length
        cds_list.append((last_start, last_length)) #optional
        try:
            new_portion = get_seq(last_chrom, seq_start, seq_end-seq_start, reference_index)
            hap_output = hap.returnphasing(last_chrom, seq_start, seq_end-1, new_portion, args.vcf)
            hap_seq_list = hap_output[0]
            for hap_seq in hap_seq_list:
                mutation_seq = ""
                wild_seq = ""
                for cds_stretch in cds_list:
                    (stretch_start, stretch_length) = cds_stretch
                    index_start = stretch_start - seq_start
                    mutation_seq += hap_seq[index_start:index_start+stretch_length]
                    wild_seq += get_seq(last_chrom, stretch_start, stretch_length, reference_index)
                kmer(mutation_posits, _size_min, _size_max,
                    seq_to_peptide(wild_seq, orf_dict[trans_id][0][0]),
                    seq_to_peptide(mutation_seq, orf_dict[trans_id][0][0])
                    )
        except:
            print "find and print kmers failure"
            raise
            return
    #if no bam: 
    else:
        wild_seq = ""
        full_length = 0
        for cds_stretch in cds_list:
            (seq_start, seq_length) = cds_stretch
            try:
                wild_seq += get_seq(last_chrom, seq_start, seq_length, reference_index)
            except:
                return
            full_length += seq_length
        cds_start = cds_list[0][0]
        #if direct == "-":
        #    mutation_seq = make_mutation_seq(wild_seq, mutation_locs, True)
        #else:
        mutation_seq = make_mutation_seq(wild_seq, mutation_locs, False)
        print(wild_seq, len(wild_seq))
        print(mutation_seq, len(mutation_seq))
        print mutation_locs
        if orf_dict[trans_id][0][0] == '-':
            reverse_strand = True
        else:
            reverse_strand = False
        kmer(mutation_posits, _size_min, _size_max,
            seq_to_peptide(wild_seq, reverse_strand=reverse_strand), 
            seq_to_peptide(mutation_seq, reverse_strand=reverse_strand)
            )

# def check_stop():
#     #Note: check slack for how to call return hapcut phasing. Edit get
#            find_seq_and_kmer, cause that is not calling it correctly rn.
#     #Needs st_ind, end_ind, the hapcut return line named "line", reference_index, get_cds stuff
#     #chrom, find_stop stuff
#     com_codon = len(line) % 3
#     #complete_codon + 1 so that you know where to have the start (chrom position) for find_stop
#     if reverse:
#         (l_side, r_side) = (com_codon+1, 0)
#     else:
#         (l_side, r_side) = (0, com_codon+1)
#     (cds_list, t1, t2) = get_cds(trans_id, [(st_ind, fake_line_number)], l_side, r_side, cds_dict, mutation_locs)
#     #NOTE: don't want the st_ind itself included
#     finish_cod = ""
#     for bound_st, bound_len in cds_list:
#         finish_cod += get_seq(chrom, bound_st, bound_len, reference_index)
#     line = finish_cod[1:] + line
#     #Now run find_stop
#     #For the wild_seq:
#     #   If forward strand: Use st_ind, find len(mutation_seq), use get_cds to find 
#     #                      that many to the right from the st_ind, use get_seq
#     #                      to construct the actual seq, print output
#     #   If reverse strand: Use end_ind, find len(mutation_seq), do opposite of above


 

def kmerize_trans(trans_lines, line_count, trans_id, trans_cds_list, direct, reference_index):
    last_chrom = "None"
    orig_seq = ""
    mutation_locs = {}
    mutation_posits = []
    mutation_line_num = 0
    bounds_set = set()
    line_num = -1
    while line_num < len(trans_lines):
        line_num += 1
        if(mutation_line_num != 0):
            line_num += mutation_line_num-1
            mutation_line_num = 0
        if(line_num >= len(trans_lines)):
            break
        line = trans_lines[line_num]
        line_count += 1
        vals = line.strip().split('\t')
        (chrom, pos, orig, alt, info) = (vals[0], int(vals[1]), vals[3], vals[4], vals[7]
            )
        tokens = info.strip().split('|')
        mutation_type = tokens[1]
        if(mutation_type != "missense_variant" and len(orig) == len(alt)): 
            continue
        if mutation_type == "missense_variant":
            rel_pos = int(tokens[13])
            pos_in_codon = (rel_pos+2)%3 #ie: ATG --> 0,1,2
        try:
            if orf_dict[trans_id][0][0] == "-" and mutation_type == "missense_variant": 
                pos_in_codon = 2-pos_in_codon
        except:
            continue
        #(cds_list, temp, bounds_set) = get_cds(trans_id, [(last_pos, line_count)], left_side, right_side, my_dict, mutation_locs)
        if last_chrom != "None":
            try:
                if direct == "+":
                    seq_end_cds_list,temp1,temp2 = get_cds(trans_id, [(last_pos, None)], 0, (3*_size_max-1)-pos_in_codon, cds_dict, {})
                    seq_end = seq_end_cds_list[-1][0] + seq_end_cds_list[-1][1] - 1
                else:
                    seq_end_cds_list,temp1,temp2 = get_cds(trans_id, [(last_pos, None)], (_size_max-1)*3+pos_in_codon, 0, cds_dict, {})
                    seq_end = seq_end_cds_list[0][0]
            except Exception as ex:
                print(type(ex))
                print "UNBOUND LOCAL ERROR 386"
                seq_end = None
        if((last_chrom != "None") and (mutation_type == "missense_variant") and seq_end != None and ((direct == "-" and pos < seq_end) or (direct == "+" and pos > seq_end))):
            print "INSIdE THE SEQ KMER THROUGHPUT"
            #(left_side,right_side) = (last_pos-st_ind,end_ind-last_pos)
            if(direct == "+"):
                left_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], last_pos)-1,
                                                        trans_cds_list)
                (left_side,right_side) = (last_pos-st_ind-left_introns, end_ind-last_pos)
            else:
                right_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], last_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                (left_side,right_side) = (last_pos-st_ind, end_ind-last_pos-right_introns)
                old_locs = mutation_locs.copy()
                mutation_locs.clear()
                for key in old_locs:
                    tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                2*bisect.bisect(trans_cds_list[::2], key+st_ind)-1,
                                                trans_cds_list)
                    mutation_locs[key-tot_introns] = old_locs[key]
            #print left_side, right_side, last_pos
            (cds_list, mutation_locs, bounds_set) = get_cds(trans_id, mutation_posits, left_side, right_side, my_dict, mutation_locs)
            if(len(cds_list) != 0):
                find_seq_and_kmer(cds_list, last_chrom, reference_index,
                                  mutation_locs, orf_dict, trans_id, mutation_posits, direct)
            (mutation_locs, mutation_posits) = (dict(), [])
        if len(orig) != len(alt):
            try:
                cds_list = cds_dict[trans_id]
            except KeyError:
                continue
            orf_list = orf_dict[trans_id]
            cds_index = 2*bisect.bisect(cds_list[::2], pos)-2
            if cds_index < 0:
                continue
            (cds_start, cds_end) = (cds_list[cds_index], cds_list[cds_index+1]) 
            orf = orf_list[cds_index//2]
            (strand, frame) = (orf[0], orf[1])
            pos_in_codon = (pos - cds_start - int(frame))%3 #check math
            if(strand != "+"):
                #print "HERE"
                #Check math.
                #print int(frame), cds_end, pos, len(alt), len(orig)
                pos_in_codon = 2-(int(frame)-(cds_end-pos+len(alt)-len(orig))%3)%3
                #pos_in_codon = (3 - (((cds_end-(pos+1)) - int(frame))%3))%3
                #print "reverse pos in codon", pos_in_codon, pos
                #pos_in_codon = 2-pos_in_codon
        if len(orig) != len(alt):
            shift = len(alt)-len(orig)
            mutation_posits.append((pos, line_count))
            if strand == "-":
                if len(mutation_locs) == 0:
                    end_ind = pos+pos_in_codon-shift+(_size_max-1)*3
                    st_ind = pos - (2-pos_in_codon)
                else:
                    old_start = st_ind
                    st_ind = pos - (2-pos_in_codon)
                    old_locs = mutation_locs.copy()
                    mutation_locs.clear()
                    for key in old_locs:
                        mutation_locs[key+old_start-st_ind] = old_locs[key]
                    #mutation_locs = dict()
                #st_ind = pos - (2-pos_in_codon)
                #@TODO don't need to set the end_ind if a missense beforehand
                #end_ind = 32 - pos_in_codon + shift + 1 + pos
                ######################################################## CHECK THIS MATH (SPECIFICALLY, THE +1)
                #st_ind = pos-(2-(pos_in_codon+abs(shift))%3)
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                query_st = st_ind
                if len(alt) > len(orig):
                    #end_ind = pos + 32 - (pos_in_codon+shift)%3
                    mutation_locs[pos-st_ind] = alt
                else:
                    #end_ind = pos + 32 - pos_in_codon + abs(shift)
                    for index in range(abs(shift)):
                        mutation_locs[pos-st_ind+1+index] = ""
                #print query_st, st_ind, end_ind, mutation_locs
                (left_side, right_side) = (pos-st_ind, end_ind-pos)
                (cds_list, mutation_locs, bounds_set) = get_cds(trans_id, mutation_posits, left_side, right_side, cds_dict, mutation_locs)
                if len(cds_list) != 0:
                    upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], cds_list[0][0])-1
                    orig_seq = ""
                    print 'cds list curr', cds_list
                    for cds_stretch in cds_list:
                        (seq_start, seq_length) = cds_stretch
                        try:
                            orig_seq += get_seq(chrom, seq_start, seq_length, reference_index)
                            print 'now', chrom, seq_start, seq_length, orig_seq
                        except:
                            print(chrom, seq_start, seq_length)
                            break
                    print 'before try', orig_seq
                    try:
                        #right_half = end_ind - pos
                        left_half, seq_st_pos = find_stop(cds_list[0][0], trans_id,
                                              line_count, cds_dict, chrom,
                                              reference_index, mutation_locs, True)
                        #print "first left", left_half
                        #print "first orig", orig_seq
                        print 'no idea where prev orig', orig_seq
                        orig_seq = left_half + orig_seq
                        print 'still no idea after orig', orig_seq
                        old_locs = mutation_locs.copy()
                        mutation_locs.clear()
                        for locate in old_locs:
                            mutation_locs[locate+len(left_half)] = old_locs[locate]
                        mutation_line_num = 0
                        new_mutation_pos = pos + 1
                        tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], seq_st_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                        while(new_mutation_pos > seq_st_pos):#end_ind - len(orig_seq)-tot_introns):
                            mutation_line_num += 1
                            if (line_num + mutation_line_num) >= len(trans_lines):
                                break
                            new_mutation = trans_lines[line_num+mutation_line_num]
                            vals = new_mutation.strip().split('\t')
                            (new_mutation_pos, orig, alt, info) = (int(vals[1]), vals[3], vals[4], vals[7]
                                )
                            if(new_mutation_pos <= seq_st_pos):#end_ind-len(orig_seq)):
                                break
                            tokens = info.strip().split('|')
                            mutation_type = tokens[1]
                            if(mutation_type != "missense_variant" and (len(orig) == len(alt))):
                                continue
                            if len(orig) != len(alt):
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], seq_st_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], new_mutation_pos)-1,
                                                        trans_cds_list)
                                print get_seq(chrom, seq_st_pos, 20, reference_index), orig_seq
                                splice_ind = new_mutation_pos-seq_st_pos-tot_introns - 1
                                print 50000000, new_mutation_pos, seq_st_pos, tot_introns, orig_seq
                                print orig_seq[splice_ind:]
                                mutation_posits.append((new_mutation_pos, line_count))
                                pos_in_codon = 2 - (pos_in_codon + (pos-new_mutation_pos+len(alt)-len(orig))%3)%3
                                #@TODO check math for pos_in_codon
                                (cds_list, t1, t2) = get_cds(trans_id, [(new_mutation_pos, line_count)], pos_in_codon, 0, cds_dict, mutation_locs)
                                shift += len(alt) - len(orig)
                                st_ind = cds_list[0][0]
                                #st_ind = new_mutation_pos - pos_in_codon
                                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                                new_left, seq_st_pos = find_stop(st_ind, trans_id, line_count, cds_dict, chrom, reference_index, mutation_locs, True)
                                prev_leng = len(orig_seq)
                                #@TODO before get_seq, run get_cds for that seq
                                #@TODO find tot_introns between seq_st_pos and end_ind, might be new_mutation_pos to end_ind
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], new_mutation_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                                #print('orig_seq is being constructed with', pos_in_codon)
                                #print('new left', new_left)
                                #print('tot_introns', tot_introns)
                                #print('get_seq', st_ind, pos_in_codon, get_seq(chrom, st_ind, pos_in_codon, reference_index))
                                #print('spliced orig', orig_seq[new_mutation_pos - (end_ind-len(orig_seq))+tot_introns-1:])
                                #print('normal prev orig', orig_seq)
                                #orig_seq = (new_left + get_seq(chrom, st_ind, pos_in_codon, reference_index) + right_half)
                                #orig_seq = (new_left + get_seq(chrom, st_ind, pos_in_codon, reference_index) +
                                #            orig_seq[new_mutation_pos - (end_ind-len(orig_seq))+tot_introns+shift:])
                                #print(len(orig_seq), 'before orig_seq len')
                                #print new_mutation_pos - (end_ind-len(orig_seq))+tot_introns-1, orig_seq[new_mutation_pos - (end_ind-len(orig_seq))+tot_introns-1:]
                                (cds_list, t1, t2) = get_cds(trans_id, [(st_ind, line_count)], 0, pos_in_codon, cds_dict, mutation_locs)
                                if len(cds_list) != 0:
                                    if len(cds_list) != 1:
                                        cds_list = cds_list[:-1] + [(cds_list[-1][0], cds_list[-1][1]-1)]
                                    else:
                                        cds_list = [(cds_list[-1][0], cds_list[-1][1]-1)]
                                mutation_codon = ""
                                for bound_start, bound_leng in cds_list:
                                    mutation_codon += get_seq(chrom, bound_start, bound_leng, reference_index)
                                print "Comparison between the two", mutation_codon, get_seq(chrom, st_ind, pos_in_codon, reference_index), pos_in_codon
                                print 'prev orig seq', orig_seq
                                print 'new left', new_left
                                print 'mutation_codon', mutation_codon
                                print 'splice', orig_seq[splice_ind:]
                                print 'mutation constructed with', cds_list
                                orig_seq = (new_left + mutation_codon + #get_seq(chrom, st_ind, pos_in_codon, reference_index) +
                                            orig_seq[splice_ind:]) #new_mutation_pos - (end_ind-len(orig_seq))+tot_introns-1:])
                                print 'by the time we are here, orig seq is', orig_seq
                                old_locs = mutation_locs.copy()
                                mutation_locs.clear()
                                for locate in old_locs:
                                    #print locate, prev_leng, len(orig_seq), locate-prev_leng+len(orig_seq)
                                    mutation_locs[locate-prev_leng+len(orig_seq)] = old_locs[locate]
                                if len(alt) > len(orig):
                                    mutation_locs[new_mutation_pos-(st_ind-len(new_left))] = alt
                                else:
                                    for delet in range(len(orig)-len(alt)):
                                        mutation_locs[new_mutation_pos+1+delet-(st_ind-len(new_left))] = ""
                            else:
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], new_mutation_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                                mutation_posits.append((new_mutation_pos, line_count))
                                print 'before mutation locs', mutation_locs
                                mutation_locs[new_mutation_pos-(end_ind-len(orig_seq))+tot_introns-1] = alt
                                print 'just added', mutation_locs
                        #@TODO set up marker for whether bam exists
                        complete_cds_list = cds_dict[trans_id]
                        lo = bisect.bisect(complete_cds_list[::2], seq_st_pos)
                        hi = bisect.bisect(complete_cds_list[1::2], end_ind) + 1
                        tupled_list = []
                        for index in range(len(complete_cds_list[lo:hi])//2):
                            tupled_list.append((complete_cds_list[2*index], complete_cds_list[2*index+1]))
                        print "BEFORE FULL_SEQ QUERY", mutation_posits, mutation_locs
                        full_seq = get_seq(chrom, seq_st_pos, end_ind-seq_st_pos, reference_index)
                        print "HAPCUT INPUT VARIABLES REVERSE < MUST INCLUDE", seq_st_pos, end_ind, tupled_list, full_seq
                        mutation_seq = make_mutation_seq(orig_seq, mutation_locs, False)
                        #wild_seq = get_seq(chrom, end_ind-len(mutation_seq)+1, len(mutation_seq), reference_index)
                        kmer(mutation_posits, _size_min, _size_max, seq_to_peptide(orig_seq, reverse_strand=True), seq_to_peptide(mutation_seq, reverse_strand=True))
                        print "Reverse Indel ", orig_seq, "\t", mutation_seq, len(orig_seq), len(mutation_seq), pos
                        print(mutation_locs, mutation_posits)
                    except Exception as ex:
                        (mutation_locs, mutation_posits, last_chrom) = (dict(), [], "None")
                        print "Reverse Failure"
                        print type(ex)
                        break
            if strand == "+":
                #pos, pos_in_codon, shift, mutation_locs, alt, orig, chrom, reference_index, line_count, exon_dict
                if(len(mutation_locs)==0):
                    st_ind = pos-(_size_max-1)*3-pos_in_codon
                    upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                if len(alt) > len(orig):
                    end_ind = pos + 2 - (pos_in_codon+shift)%3
                    #query_st = end_ind + 1
                    mutation_locs[pos-st_ind] = alt
                else:
                    end_ind = pos + 2 - pos_in_codon + abs(shift)
                    #query_st = end_ind + 1
                    for index in range(abs(shift)):
                        mutation_locs[pos-st_ind+1+index] = ""
                (left_side, right_side) = (pos-st_ind, end_ind-pos)
                (cds_list, mutation_locs, bounds_set) = get_cds(trans_id, mutation_posits, left_side, right_side, cds_dict, mutation_locs)
                if len(cds_list) != 0:
                    orig_seq = ""
                    new_start = cds_list[0][0]
                    for cds_stretch in cds_list:
                        (seq_start, seq_length) = cds_stretch
                        try:
                            orig_seq += get_seq(chrom, seq_start, seq_length, reference_index)
                        except:
                            print(chrom, seq_start, seq_length)
                            break
                    try:
                        tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                        new, seq_end_pos = find_stop(end_ind+1+tot_introns, trans_id,
                                              line_count, cds_dict, chrom,
                                              reference_index, mutation_locs, False)
                        orig_seq += new
                        mutation_line_num = new_mutation_pos = 0
                        while(new_mutation_pos < st_ind + len(orig_seq)+tot_introns):
                            mutation_line_num += 1
                            if (line_num + mutation_line_num) >= len(trans_lines):
                                break
                            new_mutation = trans_lines[line_num+mutation_line_num]
                            vals = new_mutation.strip().split('\t')
                            (new_mutation_pos, orig, alt, info) = (int(vals[1]), vals[3], vals[4], vals[7]
                                )
                            if(new_mutation_pos >= st_ind + len(orig_seq)):
                                break
                            tokens = info.strip().split('|')
                            mutation_type = tokens[1]
                            if(mutation_type != "missense_variant" and (len(orig) == len(alt))):
                                continue
                            if len(orig) != len(alt):
                                mutation_posits.append((new_mutation_pos, line_count))
                                #pos_in_codon = (new_mutation_pos-st_ind+shift)%3
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                                2*bisect.bisect(trans_cds_list[::2], new_mutation_pos)-1,
                                                                trans_cds_list)
                                if len(alt) > len(orig):
                                    mutation_locs[new_mutation_pos-st_ind-tot_introns] = alt
                                    #mutation_locs[new_mutation_pos-st_ind+shift] = alt
                                    shift += (len(alt) - len(orig))
                                    end_ind = new_mutation_pos + 2 - (pos_in_codon+shift)%3
                                    query_st = end_ind + 1
                                else:
                                    shift += (len(alt) - len(orig))
                                    end_ind = new_mutation_pos + 2 - pos_in_codon + abs(shift)
                                    query_st = end_ind + 1
                                    for index in range(abs(len(alt)-len(orig))):
                                        mutation_locs[new_mutation_pos-st_ind+1+index-tot_introns] = ""
                                pos_in_codon = (new_mutation_pos-st_ind+shift)%3
                                (cds_list, temp, bounds_set) = get_cds(trans_id, [(new_mutation_pos, line_count)], 0, 2-pos_in_codon, cds_dict, dict())
                                if len(cds_list)==0: continue
                                codon = ""
                                #print('pos in codon', pos_in_codon)
                                #print('cds list', cds_list)
                                for bound_start, bound_length in cds_list:
                                    codon += get_seq(chrom, bound_start, bound_length, reference_index)
                                final_st, final_leng = cds_list.pop()
                                #there was a +1
                                #print('before slicing and adding codon', len(orig_seq))
                                orig_seq = orig_seq[0:new_mutation_pos-st_ind-tot_introns] + codon #get_seq(chrom, new_mutation_pos, 2-pos_in_codon, reference_index)
                                new, seq_end_pos = find_stop(final_st+final_leng, trans_id, line_count, cds_dict, chrom, reference_index, mutation_locs, False)
                                orig_seq += new
                            else:
                                ################ @TODO fix line_count to be specific for mutation
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], new_mutation_pos)-1,
                                                        trans_cds_list)
                                mutation_posits.append((new_mutation_pos, line_count))
                                mutation_locs[pos-st_ind+shift-tot_introns] = alt
                        print()
                        print "HERE"
                        complete_cds_list = cds_dict[trans_id]
                        lo = bisect.bisect(complete_cds_list[::2], st_ind)
                        hi = bisect.bisect(complete_cds_list[1::2], seq_end_pos) + 1
                        tupled_list = []
                        for index in range(len(complete_cds_list[lo:hi])//2):
                            tupled_list.append((complete_cds_list[2*index], complete_cds_list[2*index+1]))
                        print "BEFORE FULL_SEQ QUERY"
                        full_seq = get_seq(chrom, st_ind, seq_end_pos-st_ind, reference_index)
                        print('ending orig_seq length is equal to ', len(orig_seq))
                        print "HAPCUT INPUT VARIABLES FORWARD < MUST INCLUDE", st_ind, seq_end_pos, tupled_list, full_seq
                        mutation_seq = make_mutation_seq(orig_seq, mutation_locs, False)
                        #print "Forward ", orig_seq, str(len(orig_seq)), str(len(mutation_seq))
                        wild_seq = get_seq(chrom, st_ind, len(mutation_seq), reference_index)
                        print(wild_seq)
                        print(mutation_seq)
                        print("len of wild_seq", len(wild_seq))
                        print("len of mutation_seq", len(mutation_seq))
                        kmer(mutation_posits, seq_to_peptide(wild_seq), seq_to_peptide(mutation_seq))
                        print "Indel ", wild_seq, "\t", mutation_seq, len(wild_seq), len(mutation_seq), pos
                        print mutation_locs
                    except KeyError:
                        pass
            (mutation_locs, mutation_posits) = (dict(), [])
            last_chrom = "None"
            continue
        if (direct == "+"):
            if last_chrom == chrom and seq_end != None and pos <= seq_end:
                #The order of that if-statement is important! Don't change it!
                end_ind = pos+(3*_size_max-1)-pos_in_codon
            else:
                st_ind = pos-(_size_max-1)*3-pos_in_codon
                #(cds_list, temp, temp2) = get_cds(trans_id, [(pos, line_count)], pos-30-pos_in_codon, 0, cds_dict, mutation_locs)
                #st_ind = cds_list[0][0] + cds_list[0][1]
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                end_ind = pos+(3*_size_max-1)-pos_in_codon
            mutation_locs[(pos-st_ind)] = alt
        else:
            if last_chrom == chrom and seq_end != None and pos >= seq_end:
                #The order of that if-statement is important! Don't change it!
                old_start = st_ind
                st_ind = pos-(_size_max-1)*3-pos_in_codon
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                old_locs = mutation_locs.copy()
                mutation_locs.clear()
                for key in old_locs:
                    mutation_locs[key+old_start-st_ind] = old_locs[key]
            else:
                #@TODO check
                st_ind = pos-(_size_max-1)*3-pos_in_codon
                #(cds_list, temp, temp2) = get_cds(trans_id, [(pos, line_count)], 0, pos+32-pos_in_codon, cds_dict, mutation_locs)
                #end_ind = cds_list[-1][0] + cds_list[-1][1] - 1
                end_ind = pos+(3*_size_max-1)-pos_in_codon
                print "696 after calculations were made here is here", pos, st_ind, end_ind
            mutation_locs[pos-st_ind] = alt
        mutation_posits.append((pos, line_count))
        (last_pos,last_chrom) = (pos, chrom)
        print "End of loop", mutation_posits
    print "outside for loop", mutation_posits
    try:
        if(direct == "+"):
            left_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                    2*bisect.bisect(trans_cds_list[::2], last_pos)-1,
                                                    trans_cds_list)
            (left_side,right_side) = (last_pos-st_ind-left_introns, end_ind-last_pos)
        else:
            right_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], last_pos)-1,
                                                    2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                    trans_cds_list)
            (left_side,right_side) = (last_pos-st_ind, end_ind-last_pos-right_introns)
            old_locs = mutation_locs.copy()
            mutation_locs.clear()
            for key in old_locs:
                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                            2*bisect.bisect(trans_cds_list[::2], key+st_ind)-1,
                                            trans_cds_list)
                print "here", key, tot_introns
                mutation_locs[key-tot_introns] = old_locs[key]
    except:
        return
    (cds_list,mutation_locs, bounds_set) = get_cds(trans_id, mutation_posits, left_side, right_side, my_dict, mutation_locs)
    print "cds_list", len(cds_list), "st_ind", st_ind, "end_ind", end_ind, "mutation posits", mutation_posits
    #print('cds', cds_list)
    #print(mutation_posits, left_side, right_side)
    #print("Final mutation_locs: ", str(mutation_locs), str(mutation_posits))
    if(len(cds_list) != 0):
        find_seq_and_kmer(cds_list, last_chrom, reference_index, mutation_locs,
                          orf_dict, trans_id, mutation_posits, direct)


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
