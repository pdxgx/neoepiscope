import bisect
import argparse
import bowtie_index
import sys
import math
import string
import copy
import pickle
#import Hapcut2interpreter as hap
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
args = parser.parse_args()
ref_ind = bowtie_index.BowtieIndexReference(args.bowtie_index)

if args.bam:
    try:
        import Hapcut2interpreter as hap
    except:
        raise RuntimeError('Hapcut2interpreter was not run')

my_dicts = pickle.load ( open (args.dicts, "rb"))
(cds_dict, orf_dict, exon_dict, exon_orf_dict) = (my_dicts[0], my_dicts[1], my_dicts[2], my_dicts[3])

codon_table = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"Stop", "TAG":"Stop",
    "TGT":"C", "TGC":"C", "TGA":"Stop", "TGG":"W",
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
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def turn_to_aa(nucleotide_string, strand="+"):
    ''' References codon_table to turn codons to amino acids.
        nucleotide_string: (String) Entire nucleotide sequence made from A,T,C,G
        strand: (String) Denotes forward or reverse strand
        Return value: (String) Amino acid sequence 
    '''
    aa_string = ""
    num_aa = 0
    if strand == "-":
        translation_table = string.maketrans("ATCG", "TAGC")
        nucleotide_string = nucleotide_string.translate(translation_table)[::-1]
    for aa in range(len(nucleotide_string)//3):
        num_aa += 1
        try:
            codon = codon_table[nucleotide_string[3*aa:3*aa+3]]
        except KeyError:
            print >>sys.stderr, (
                        'Could not translate nucleotide string "{}".'
                    ).format(nucleotide_string)
            return False
        if codon == "Stop":
            aa_string += (len(nucleotide_string)//3 - num_aa + 1)*'X'
            break
        else:
            aa_string += codon
    return aa_string

def my_print_function(kmer_list, mute_posits):
    ''' Prints out epitopes in 3 columns: Wild Type, Mutant Type, and
            Position on chromosome
        kmer_list: (List) Tuples containing Wild Types and Mutant Types for
            printing
        mute_posits: (List) Positions of mutations on chromosome
        Return value: (String) Amino acid sequence 
    '''
    if len(kmer_list)==0: return None
    for wtmtPair in kmer_list:
        wt,mt = wtmtPair
        print(wt + "\t" + mt + "\t" + str(mute_posits))
    return None


def kmer(mute_posits, normal_aa, mutated_aa = ""):
    ''' Generates 8-11' kmers from amino acid sequences.
        mute_posits: (List) Positions of mutations on chromosome
        normal_aa: (String) Amino acid sequence of Wild Type
        mutated_aa: (String) Amino acid sequence of Mutant Type
        Return value: (List) List of tuples pairing mismatching kmers for a
            given 8-11' epitope.
    '''
    if (len(mutated_aa) == 0):
        mutated_aa = normal_aa
    kmer_list = list()
    #Loop through window sizes
    for ksize in range(8, 12):
        for startIndex in range(len(mutated_aa)-ksize+1):
            kmer_list.append((normal_aa[startIndex:startIndex+ksize], mutated_aa[startIndex:startIndex+ksize]))
    final_list = list()
    for WT,MT in kmer_list:
        if (WT != MT and 'X' not in MT):
            final_list.append((WT, MT))
    my_print_function(final_list, mute_posits)
    return final_list

def get_cds(transcript_id, mutation_pos_list, seq_length_left, 
              seq_length_right, ordered_cds_dict, mute_dict):
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
        sequence necessary for 8-11' peptide kmerization based on the position 
        of a mutation within a chromosome.
    '''
    bounds_set = set()
    if transcript_id not in ordered_cds_dict:
        print('failed on 1')
        print(transcript_id)
        print(ordered_cds_dict)
        return [], mute_dict, bounds_set
    pos_in_codon = 2 - (seq_length_right%3)
    cds_list = ordered_cds_dict[transcript_id]
    mutation_pos = -1
    #Don't want to check rightmost since seq. queries based off of it.
    if len(mutation_pos_list) >= 2:
        removal_list = []
        shift = mutation_pos_list[0][0] - min(mute_dict)
        key_list = list(mute_dict.keys())
        #Remove all mutations outside of cds boundaries.
        for index in range(len(mutation_pos_list)):
            lower_cds_index = 2*bisect.bisect(cds_list[::2], mutation_pos_list[index][0])-2
            upper_cds_index = lower_cds_index+1
            if(lower_cds_index < 0 or 
               cds_list[upper_cds_index] < mutation_pos_list[index][0]):
                #Delete at the current index
                try:
                    if mute_dict[key_list[index]] == '':
                        try:
                            count = 1
                            while mute_dict[key_list[index+count]] == '':
                                count += 1
                                removal_list.append(index+count)
                                del mute_dict[key_list[index+count]]
                        except IndexError:
                            pass
                    del mute_dict[key_list[index]]
                    removal_list.append(index)
                except IndexError:
                    continue
        for index in range(len(removal_list)-1, -1, -1):
            #print("made edits to mute pos list")
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
            return [], mute_dict, bounds_set
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
                seq_length_right = 30 + new_pos_in_codon
                seq_length_left -= (mutation_pos_list[-1][0] 
                                    - mutation_pos_list[index][0])
            break
    if(mutation_pos == -1):
        print('failed on 3')
        return [], mute_dict, bounds_set
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
    return nucleotide_index_list, mute_dict, bounds_set

def get_introns(start_index, stop_index, cds_list):
    total_introns = 0
    index = start_index
    while index < stop_index:
        total_introns += cds_list[index+1]-cds_list[index]-1
        index += 2
    return total_introns

def find_stop(query_st, trans_id, line_count, cds_dict, chrom, ref_ind, mute_locs, reverse):
    ''' Queries get_cds() and Bowtie to find a stop codon in the case of a phase
            shift (indel)
        query_st: (int)
        trans_id: ()
        line_count: ()
        cds_dict: ()
        chrom: (int) The chromosome that the mutation is located on
        ref_ind: ()
        mute_locs: (Dict) Dictionary mapping sequence position to 
            actual mutation
        reverse: (Bool) True if mutation represented by reverse strand
        Return value:  
    '''
    until_stop = ""
    start = query_st
    stop_found = False
    (l_query, r_query) = (0, 33)
    if reverse:
        (l_query, r_query) = (r_query, l_query)
    while(stop_found == False):
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(start,line_count)], l_query, r_query, cds_dict, mute_locs)
        extra_cods = ""
        for bound_start, bound_stretch in exon_list:
            extra_cods += get_seq(chrom, bound_start, bound_stretch, ref_ind)
        if reverse:
        #    extra_cods = extra_cods[:34]
            start = exon_list[0][0]
        #    print "extra_cods", start, query_st, l_query, extra_cods
        else:
            start = exon_list[-1][0] + exon_list[-1][1] - 1
        count = 0
        while(count<33):
            if reverse:
                new_codon = extra_cods[30-count:33-count]
                count += 3
                amino_acid = turn_to_aa(new_codon, "-")
            else:
                new_codon = extra_cods[count: count+3]
                count += 3
                amino_acid = turn_to_aa(new_codon)
            if(amino_acid==""):
                stop_found = True
                break
            if reverse:
                until_stop = new_codon + until_stop
            else:
                until_stop += new_codon
    if(reverse):
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(start,line_count)], count, 0, cds_dict, mute_locs)
        return until_stop, exon_list[0][0]
    else:
        (exon_list, temp_out, bounds_set) = get_cds(trans_id, [(start,line_count)], 0, count, cds_dict, mute_locs)
        return until_stop, exon_list[-1][0]+exon_list[-1][1]-1
    #return until_stop

def get_seq(chrom, start, splice_length, ref_ind):
    chr_name = "chr" + chrom #proper
    start -= 1 #adjust for 0-based bowtie queries
    try:
        seq = ref_ind.get_stretch(chr_name, start, splice_length)
    except KeyError:
        return False
    return seq

def make_mute_seq(orig_seq, mute_locs, reverse):
    mute_seq = ""
    if(reverse):
        orig_seq = orig_seq[::1]
        for ind in range(len(orig_seq)):
            if ind in mute_locs:
                mute_seq = mute_locs[ind]+mute_seq
            else:
                mute_seq = orig_seq[ind] + mute_seq 
        return mute_seq[::-1]
    else:
        for ind in range(len(orig_seq)):
            if ind in mute_locs:
                mute_seq += mute_locs[ind]
            else:
                mute_seq += orig_seq[ind]
        return mute_seq

def find_seq_and_kmer(cds_list, last_chrom, ref_ind, mute_locs,
                      orf_dict, trans_id, mute_posits, direct):
    #if bam exists:
    if args.bam:
        seq_start = cds_list[0][0]
        (last_start, last_length) = cds_list.pop()
        seq_end = last_start+last_length
        cds_list.append((last_start, last_length)) #optional
        try:
            new_portion = get_seq(last_chrom, seq_start, seq_end-seq_start, ref_ind)
            hap_output = hap.returnphasing(last_chrom, seq_start, seq_end-1, new_portion, args.vcf)
            hap_seq_list = hap_output[0]
            for hap_seq in hap_seq_list:
                mute_seq = ""
                wild_seq = ""
                for cds_stretch in cds_list:
                    (stretch_start, stretch_length) = cds_stretch
                    index_start = stretch_start - seq_start
                    mute_seq += hap_seq[index_start:index_start+stretch_length]
                    wild_seq += get_seq(last_chrom, stretch_start, stretch_length, ref_ind)
                kmer(mute_posits,
                    turn_to_aa(wild_seq, orf_dict[trans_id][0][0]),
                    turn_to_aa(mute_seq, orf_dict[trans_id][0][0])
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
                wild_seq += get_seq(last_chrom, seq_start, seq_length, ref_ind)
            except:
                return
            full_length += seq_length
        cds_start = cds_list[0][0]
        if direct == "-":
            mute_seq = make_mute_seq(wild_seq, mute_locs, True)
        else:
            mute_seq = make_mute_seq(wild_seq, mute_locs, False)
        print(wild_seq, len(wild_seq))
        print(mute_seq, len(mute_seq))
        kmer(mute_posits,
            turn_to_aa(wild_seq, orf_dict[trans_id][0][0]), 
            turn_to_aa(mute_seq, orf_dict[trans_id][0][0])
            )
    

def kmerize_trans(trans_lines, line_count, trans_id, trans_cds_list, direct):
    last_chrom = "None"
    orig_seq = ""
    mute_locs = {}
    mute_posits = []
    mute_line_num = 0
    bounds_set = set()
    line_num = -1
    while line_num < len(trans_lines):
        line_num += 1
        if(mute_line_num != 0):
            line_num += mute_line_num-1
            mute_line_num = 0
        if(line_num >= len(trans_lines)):
            break
        line = trans_lines[line_num]
        line_count += 1
        vals = line.strip().split('\t')
        (chrom, pos, orig, alt, info) = (vals[0], int(vals[1]), vals[3], vals[4], vals[7]
            )
        tokens = info.strip().split('|')
        mute_type = tokens[1]
        if(mute_type != "missense_variant" and len(orig) == len(alt)): 
            continue
        if mute_type == "missense_variant":
            rel_pos = int(tokens[13])
            pos_in_codon = (rel_pos+2)%3 #ie: ATG --> 0,1,2
        try:
            if orf_dict[trans_id][0][0] == "-" and mute_type == "missense_variant": 
                pos_in_codon = 2-pos_in_codon
        except:
            continue
        #(cds_list, temp, bounds_set) = get_cds(trans_id, [(last_pos, line_count)], left_side, right_side, my_dict, mute_locs)
        if last_chrom != "None":
            try:
                if direct == "+":
                    seq_end_cds_list,temp1,temp2 = get_cds(trans_id, [(last_pos, None)], 0, 32-pos_in_codon, cds_dict, {})
                    seq_end = seq_end_cds_list[-1][0] + seq_end_cds_list[-1][1] - 1
                else:
                    seq_end_cds_list,temp1,temp2 = get_cds(trans_id, [(last_pos, None)], 30+pos_in_codon, 0, cds_dict, {})
                    seq_end = seq_end_cds_list[0][0]
            except Exception as ex:
                print(type(ex))
                print "UNBOUND LOCAL ERROR 386"
                seq_end = None
        if((last_chrom != "None") and (mute_type == "missense_variant") and seq_end != None and ((direct == "-" and pos < seq_end) or (direct == "+" and pos > seq_end))):
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
            print left_side, right_side, last_pos
            (cds_list, mute_locs, bounds_set) = get_cds(trans_id, mute_posits, left_side, right_side, my_dict, mute_locs)
            if(len(cds_list) != 0):
                find_seq_and_kmer(cds_list, last_chrom, ref_ind,
                                  mute_locs, orf_dict, trans_id, mute_posits, direct)
            (mute_locs, mute_posits) = (dict(), [])
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
                print "HERE"
                #Check math.
                print int(frame), cds_end, pos, len(alt), len(orig)
                pos_in_codon = 2-(int(frame)-(cds_end-pos+len(alt)-len(orig))%3)%3
                #pos_in_codon = (3 - (((cds_end-(pos+1)) - int(frame))%3))%3
                print "reverse pos in codon", pos_in_codon, pos
                #pos_in_codon = 2-pos_in_codon
        if len(orig) != len(alt):
            shift = len(alt)-len(orig)
            mute_posits.append((pos, line_count))
            if strand == "-":
                if len(mute_locs) == 0:
                    end_ind = pos+pos_in_codon-shift+30
                st_ind = pos - (2-pos_in_codon)
                mute_locs = dict()
                #@TODO don't need to set the end_ind if a missense beforehand
                #end_ind = 32 - pos_in_codon + shift + 1 + pos
                print "END INDEX", end_ind, shift, pos
                ######################################################## CHECK THIS MATH (SPECIFICALLY, THE +1)
                #st_ind = pos-(2-(pos_in_codon+abs(shift))%3)
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                query_st = st_ind
                if len(alt) > len(orig):
                    #end_ind = pos + 32 - (pos_in_codon+shift)%3
                    mute_locs[pos-st_ind] = alt
                    print pos-st_ind, 'first locate'
                else:
                    #end_ind = pos + 32 - pos_in_codon + abs(shift)
                    for index in range(abs(shift)):
                        mute_locs[pos-st_ind+1+index] = ""
                        print pos-st_ind+1+index, 'other first locate'
                #print query_st, st_ind, end_ind, mute_locs
                print "reverse", str(end_ind-st_ind), str(shift), str(pos_in_codon)
                (left_side, right_side) = (pos-st_ind, end_ind-pos)
                (cds_list, mute_locs, bounds_set) = get_cds(trans_id, mute_posits, left_side, right_side, cds_dict, mute_locs)
                print 465, mute_posits
                if len(cds_list) != 0:
                    upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], cds_list[0][0])-1
                    orig_seq = ""
                    for cds_stretch in cds_list:
                        (seq_start, seq_length) = cds_stretch
                        try:
                            orig_seq += get_seq(chrom, seq_start, seq_length, ref_ind)
                        except:
                            print(chrom, seq_start, seq_length)
                            break
                    try:
                        #right_half = end_ind - pos
                        left_half, seq_st_pos = find_stop(cds_list[0][0], trans_id,
                                              line_count, cds_dict, chrom,
                                              ref_ind, mute_locs, True)
                        print "seq_st_pos", seq_st_pos, shift
                        #print "first left", left_half
                        #print "first orig", orig_seq
                        orig_seq = left_half + orig_seq
                        print 515, len(left_half), len(orig_seq)
                        old_locs = mute_locs.copy()
                        mute_locs.clear()
                        for locate in old_locs:
                            mute_locs[locate+len(left_half)] = old_locs[locate]
                            print locate, 'higher locate'
                        mute_line_num = 0
                        new_mute_pos = pos + 1
                        tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], seq_st_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                        while(new_mute_pos > seq_st_pos):#end_ind - len(orig_seq)-tot_introns):
                            print 498
                            mute_line_num += 1
                            if (line_num + mute_line_num) >= len(trans_lines):
                                break
                            new_mute = trans_lines[line_num+mute_line_num]
                            vals = new_mute.strip().split('\t')
                            (new_mute_pos, orig, alt, info) = (int(vals[1]), vals[3], vals[4], vals[7]
                                )
                            print new_mute_pos, orig, alt
                            #print new_mute_pos, end_ind-len(orig_seq)
                            if(new_mute_pos <= seq_st_pos):#end_ind-len(orig_seq)):
                                break
                            tokens = info.strip().split('|')
                            mute_type = tokens[1]
                            if(mute_type != "missense_variant" and (len(orig) == len(alt))):
                                continue
                            if len(orig) != len(alt):
                                mute_posits.append((new_mute_pos, line_count))
                                pos_in_codon = 2 - (pos_in_codon + (pos-new_mute_pos+len(alt)-len(orig))%3)%3
                                #@TODO check math for pos_in_codon
                                shift += len(alt) - len(orig)
                                st_ind = new_mute_pos - pos_in_codon
                                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                                new_left, seq_st_pos = find_stop(st_ind, trans_id, line_count, cds_dict, chrom, ref_ind, mute_locs, True)
                                prev_leng = len(orig_seq)
                                #@TODO before get_seq, run get_cds for that seq
                                #@TODO find tot_introns between seq_st_pos and end_ind, might be new_mute_pos to end_ind
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], new_mute_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                                #print('orig_seq is being constructed with', pos_in_codon)
                                print('new left', new_left)
                                print('tot_introns', tot_introns)
                                #print('get_seq', st_ind, pos_in_codon, get_seq(chrom, st_ind, pos_in_codon, ref_ind))
                                #print('spliced orig', orig_seq[new_mute_pos - (end_ind-len(orig_seq))+tot_introns-1:])
                                #print('normal prev orig', orig_seq)
                                #orig_seq = (new_left + get_seq(chrom, st_ind, pos_in_codon, ref_ind) + right_half)
                                #orig_seq = (new_left + get_seq(chrom, st_ind, pos_in_codon, ref_ind) +
                                #            orig_seq[new_mute_pos - (end_ind-len(orig_seq))+tot_introns+shift:])
                                print(len(orig_seq), 'before orig_seq len')
                                new_var = len(orig_seq)
                                orig_seq = (new_left + get_seq(chrom, st_ind, pos_in_codon, ref_ind) +
                                            orig_seq[new_mute_pos - (end_ind-len(orig_seq))+tot_introns-1:])
                                print(len(orig_seq), 'after orig_seq len')
                                print('new_mute_pos, end_ind, old len(orig_seq), tot_introns', new_mute_pos, end_ind, new_var, tot_introns)
                                print(end_ind-new_var, 'end_ind-new_var')
                                print(new_mute_pos-(end_ind-new_var), 'new_mute_pos-above')
                                print(new_mute_pos - (end_ind-new_var)+tot_introns, 'final. -1 needs to be added')
                                old_locs = mute_locs.copy()
                                mute_locs.clear()
                                print 569, len(orig_seq), old_locs
                                for locate in old_locs:
                                    #print locate, prev_leng, len(orig_seq), locate-prev_leng+len(orig_seq)
                                    mute_locs[locate-prev_leng+len(orig_seq)] = old_locs[locate]
                                    print locate, "locate"
                                    print locate-prev_leng+len(orig_seq), "new"
                                if len(alt) > len(orig):
                                    mute_locs[new_mute_pos-(st_ind-len(new_left))] = alt
                                else:
                                    for delet in range(len(orig)-len(alt)):
                                        mute_locs[new_mute_pos+1+delet-(st_ind-len(new_left))] = ""
                            else:
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], new_mute_pos)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                                print 555555555, new_mute_pos, end_ind, len(orig_seq), tot_introns
                                print new_mute_pos-(end_ind-len(orig_seq))+tot_introns
                                mute_posits.append((new_mute_pos, line_count))
                                mute_locs[new_mute_pos-(end_ind-len(orig_seq))+tot_introns-1] = alt
                                print new_mute_pos-(end_ind-len(orig_seq))+tot_introns-1, 'prob on this one'
                                print new_mute_pos, end_ind, len(orig_seq), tot_introns
                        #@TODO set up marker for whether bam exists
                        complete_cds_list = cds_dict[trans_id]
                        lo = bisect.bisect(complete_cds_list[::2], seq_st_pos)
                        hi = bisect.bisect(complete_cds_list[1::2], end_ind) + 1
                        tupled_list = []
                        for index in range(len(complete_cds_list[lo:hi])//2):
                            tupled_list.append((complete_cds_list[2*index], complete_cds_list[2*index+1]))

                        print "BEFORE FULL_SEQ QUERY", mute_posits, mute_locs
                        full_seq = get_seq(chrom, seq_st_pos, end_ind-seq_st_pos, ref_ind)
                        print "HAPCUT INPUT VARIABLES REVERSE < MUST INCLUDE", seq_st_pos, end_ind, tupled_list, full_seq
                        mute_seq = make_mute_seq(orig_seq, mute_locs, False)
                        #wild_seq = get_seq(chrom, end_ind-len(mute_seq)+1, len(mute_seq), ref_ind)
                        kmer(mute_posits, turn_to_aa(orig_seq, "-"), turn_to_aa(mute_seq, "-"))
                        print "Reverse Indel ", orig_seq, "\t", mute_seq, len(orig_seq), len(mute_seq), pos
                        print(mute_locs, mute_posits)
                    except Exception as ex:
                        (mute_locs, mute_posits, last_chrom) = (dict(), [], "None")
                        print "Reverse Failure"
                        print type(ex)
                        break
            if strand == "+":
                #pos, pos_in_codon, shift, mute_locs, alt, orig, chrom, ref_ind, line_count, exon_dict
                if(len(mute_locs)==0):
                    st_ind = pos-30-pos_in_codon
                    upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                if len(alt) > len(orig):
                    end_ind = pos + 2 - (pos_in_codon+shift)%3
                    #query_st = end_ind + 1
                    mute_locs[pos-st_ind] = alt
                else:
                    end_ind = pos + 2 - pos_in_codon + abs(shift)
                    #query_st = end_ind + 1
                    for index in range(abs(shift)):
                        mute_locs[pos-st_ind+1+index] = ""
                (left_side, right_side) = (pos-st_ind, end_ind-pos)
                (cds_list, mute_locs, bounds_set) = get_cds(trans_id, mute_posits, left_side, right_side, cds_dict, mute_locs)
                if len(cds_list) != 0:
                    orig_seq = ""
                    new_start = cds_list[0][0]
                    for cds_stretch in cds_list:
                        (seq_start, seq_length) = cds_stretch
                        try:
                            orig_seq += get_seq(chrom, seq_start, seq_length, ref_ind)
                        except:
                            print(chrom, seq_start, seq_length)
                            break
                    try:
                        tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                                        trans_cds_list)
                        new, seq_end_pos = find_stop(end_ind+1+tot_introns, trans_id,
                                              line_count, cds_dict, chrom,
                                              ref_ind, mute_locs, False)
                        orig_seq += new
                        mute_line_num = new_mute_pos = 0
                        while(new_mute_pos < st_ind + len(orig_seq)+tot_introns):
                            mute_line_num += 1
                            if (line_num + mute_line_num) >= len(trans_lines):
                                break
                            new_mute = trans_lines[line_num+mute_line_num]
                            vals = new_mute.strip().split('\t')
                            (new_mute_pos, orig, alt, info) = (int(vals[1]), vals[3], vals[4], vals[7]
                                )
                            if(new_mute_pos >= st_ind + len(orig_seq)):
                                break
                            tokens = info.strip().split('|')
                            mute_type = tokens[1]
                            if(mute_type != "missense_variant" and (len(orig) == len(alt))):
                                continue
                            if len(orig) != len(alt):
                                mute_posits.append((new_mute_pos, line_count))
                                #pos_in_codon = (new_mute_pos-st_ind+shift)%3
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                                2*bisect.bisect(trans_cds_list[::2], new_mute_pos)-1,
                                                                trans_cds_list)
                                if len(alt) > len(orig):
                                    mute_locs[new_mute_pos-st_ind-tot_introns] = alt
                                    #mute_locs[new_mute_pos-st_ind+shift] = alt
                                    shift += (len(alt) - len(orig))
                                    end_ind = new_mute_pos + 2 - (pos_in_codon+shift)%3
                                    query_st = end_ind + 1
                                else:
                                    shift += (len(alt) - len(orig))
                                    end_ind = new_mute_pos + 2 - pos_in_codon + abs(shift)
                                    query_st = end_ind + 1
                                    for index in range(abs(len(alt)-len(orig))):
                                        mute_locs[new_mute_pos-st_ind+1+index-tot_introns] = ""
                                pos_in_codon = (new_mute_pos-st_ind+shift)%3
                                (cds_list, temp, bounds_set) = get_cds(trans_id, [(new_mute_pos, line_count)], 0, 2-pos_in_codon, cds_dict, dict())
                                if len(cds_list)==0: continue
                                codon = ""
                                #print('pos in codon', pos_in_codon)
                                #print('cds list', cds_list)
                                for bound_start, bound_length in cds_list:
                                    codon += get_seq(chrom, bound_start, bound_length, ref_ind)
                                final_st, final_leng = cds_list.pop()
                                #there was a +1
                                #print('before slicing and adding codon', len(orig_seq))
                                orig_seq = orig_seq[0:new_mute_pos-st_ind-tot_introns] + codon #get_seq(chrom, new_mute_pos, 2-pos_in_codon, ref_ind)
                                #print('after slicing and adding codon', len(orig_seq))
                                #print('len of slice', new_mute_pos-st_ind-tot_introns)
                                #print('codon elngth', len(codon))
                                #new = find_stop(st_ind+len(orig_seq)+tot_introns, trans_id, line_count, cds_dict, chrom, ref_ind, mute_locs, False)
                                new, seq_end_pos = find_stop(final_st+final_leng, trans_id, line_count, cds_dict, chrom, ref_ind, mute_locs, False)
                                #tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                #                                2*bisect.bisect(trans_cds_list[::2], seq_end_pos)-1,
                                #                                trans_cds_list)
                                orig_seq += new
                                #print('length of new (from find_stop() function)', len(new))
                                #print "final part", shift, len(orig_seq)
                            else:
                                ################ @TODO fix line_count to be specific for mute
                                tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
                                                        2*bisect.bisect(trans_cds_list[::2], new_mute_pos)-1,
                                                        trans_cds_list)
                                mute_posits.append((new_mute_pos, line_count))
                                mute_locs[pos-st_ind+shift-tot_introns] = alt
                        print()
                        print "HERE"
                        complete_cds_list = cds_dict[trans_id]
                        lo = bisect.bisect(complete_cds_list[::2], st_ind)
                        hi = bisect.bisect(complete_cds_list[1::2], seq_end_pos) + 1
                        tupled_list = []
                        for index in range(len(complete_cds_list[lo:hi])//2):
                            tupled_list.append((complete_cds_list[2*index], complete_cds_list[2*index+1]))
                        print "BEFORE FULL_SEQ QUERY"
                        full_seq = get_seq(chrom, st_ind, seq_end_pos-st_ind, ref_ind)
                        print('ending orig_seq length is equal to ', len(orig_seq))
                        print "HAPCUT INPUT VARIABLES FORWARD < MUST INCLUDE", st_ind, seq_end_pos, tupled_list, full_seq
                        mute_seq = make_mute_seq(orig_seq, mute_locs, False)
                        #print "Forward ", orig_seq, str(len(orig_seq)), str(len(mute_seq))
                        wild_seq = get_seq(chrom, st_ind, len(mute_seq), ref_ind)
                        print(wild_seq)
                        print(mute_seq)
                        print("len of wild_seq", len(wild_seq))
                        print("len of mute_seq", len(mute_seq))
                        kmer(mute_posits, turn_to_aa(wild_seq, "+"), turn_to_aa(mute_seq, "+"))
                        print "Indel ", wild_seq, "\t", mute_seq, len(wild_seq), len(mute_seq), pos
                        print mute_locs
                    except KeyError:
                        pass
            (mute_locs, mute_posits) = (dict(), [])
            last_chrom = "None"
            continue
        if (direct == "+"):
            if last_chrom == chrom and seq_end != None and pos <= seq_end:
                #The order of that if-statement is important! Don't change it!
                end_ind = pos+32-pos_in_codon
            else:
                st_ind = pos-30-pos_in_codon
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                end_ind = pos+32-pos_in_codon
            mute_locs[(pos-st_ind)] = alt
        else:
            if last_chrom == chrom and seq_end != None and pos >= seq_end:
                #The order of that if-statement is important! Don't change it!
                st_ind = pos-30-pos_in_codon
                upper_st_ind_in_cds = 2*bisect.bisect(trans_cds_list[::2], st_ind)-1
                #end_ind = pos+32-pos_in_codon
            else:
                #@TODO check
                st_ind = pos-30-pos_in_codon
                end_ind = pos+32-pos_in_codon
                print "696 after calculations were made here is here", pos, st_ind, end_ind
            #tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], st_ind)-1,
            #                        2*bisect.bisect(trans_cds_list[::2], pos)-1,
            #                        trans_cds_list)
            tot_introns = get_introns(2*bisect.bisect(trans_cds_list[::2], pos)-1,
                                   2*bisect.bisect(trans_cds_list[::2], end_ind)-1,
                                   trans_cds_list)
            mute_locs[(end_ind-pos-tot_introns-1)] = alt
            #mute_locs[(pos-st_ind+tot_introns)] = alt
        #mute_locs[(pos-st_ind)] = alt
        mute_posits.append((pos, line_count))
        (last_pos,last_chrom) = (pos, chrom)
        print "End of loop", mute_posits
    print "outside for loop", mute_posits
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
    except:
        return
    (cds_list,mute_locs, bounds_set) = get_cds(trans_id, mute_posits, left_side, right_side, my_dict, mute_locs)
    print "cds_list", len(cds_list), "st_ind", st_ind, "end_ind", end_ind, "mute posits", mute_posits
    #print('cds', cds_list)
    #print(mute_posits, left_side, right_side)
    #print("Final mute_locs: ", str(mute_locs), str(mute_posits))
    if(len(cds_list) != 0):
        find_seq_and_kmer(cds_list, last_chrom, ref_ind, mute_locs,
                          orf_dict, trans_id, mute_posits, direct)




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
        #print(get_seq("3", 69990, 20, ref_ind))
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