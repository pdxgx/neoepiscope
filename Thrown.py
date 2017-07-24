import bisect
import argparse
import bowtie_index
import sys
import string
import pickle
import Hapcut2interpreter as hap

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
    aa_string = ""
    if strand == "-":
        translation_table = string.maketrans("ATCG", "TAGC")
        nucleotide_string = nucleotide_string.translate(translation_table)[::-1]
    for aa in range(len(nucleotide_string)//3):
        try:
            codon = codon_table[nucleotide_string[3*aa:3*aa+3]]
        except KeyError:
            # print >>sys.stderr, (
            #             'Could not translate nucleotide string "{}".'
            #         ).format(nucleotide_string)
            return ""
        if (codon == "Stop"):
            break
        else:
            aa_string += codon
    return aa_string

def my_print_function(kmer_list, mute_posits):
    if len(kmer_list)==0: return None
    for wtmtPair in kmer_list:
        wt,mt = wtmtPair
        print(wt + "\t" + mt + "\t" + str(mute_posits))
    return None


def kmer(mute_posits, normal_aa, mutated_aa = ""):
    if (len(mutated_aa) == 0):
        mutated_aa = normal_aa
    kmer_list = list()
    #Loop through window sizes
    for ksize in range(8, 12):
        for startIndex in range(len(mutated_aa)-ksize):
            kmer_list.append((normal_aa[startIndex:startIndex+ksize], mutated_aa[startIndex:startIndex+ksize]))
    final_list = list()
    for WT,MT in kmer_list:
        if (WT != MT):
            final_list.append((WT, MT))
    my_print_function(final_list, mute_posits)
    return final_list

def get_cds(transcript_id, mutation_pos_list, seq_length_left, 
              seq_length_right, cds_dict, mute_dict):
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
    ordered_cds_dict = cds_dict
    bounds_set = set()
    if transcript_id not in ordered_cds_dict:
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
                    del mute_dict[key_list[index]]
                    removal_list.append(index)
                except KeyError:
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


def get_seq(chrom, start, splice_length, ref_ind):
    chr_name = "chr" + chrom
    start -= 1 #adjust for 0-based bowtie queries
    try:
        seq = ref_ind.get_stretch(chr_name, start, splice_length)
    except KeyError:
        return False
    return seq

def make_mute_seq(orig_seq, mute_locs):
    mute_seq = ""
    for ind in range(len(orig_seq)):
        if ind in mute_locs:
            mute_seq += mute_locs[ind]
        else:
            mute_seq += orig_seq[ind]
    return mute_seq

def find_seq_and_kmer(cds_list, last_chrom, ref_ind, mute_locs,
                      orf_dict, trans_id, mute_posits):
    hap_seq_list = []
    #if bam exists:
    seq_start = cds_list[0][0]
    (last_start, last_length) = cds_list.pop()
    seq_end = last_start+last_length
    cds_list.append((last_start, last_length)) #optional
    try:
        new_portion = get_seq(last_chrom, seq_start, seq_end-seq_start, ref_ind)
        hap_output = hap.returnphasing(last_chrom, seq_start, seq_end-1, new_portion, args.vcf)
        if((len(hap_output) == 2) or (len(hap_output) == 5)):
            hap_seq_list.append(hap_output[0])
        else:
            hap_seq_list.append(hap_output[0])
            hap_seq_list.append(hap_output[2])
        #print(hap_seq_list)
        for hap_seq in hap_seq_list:
            mute_seq = ""
            wild_seq = ""
            for cds_stretch in cds_list:
                (stretch_start, stretch_length) = cds_stretch
                index_start = stretch_start - seq_start
                mute_seq += hap_seq[index_start:index_start+stretch_length]
                wild_seq += get_seq(last_chrom, stretch_start, stretch_length, ref_ind)
            #mute_seq = make_mute_seq(wild_seq, mute_locs)
            kmer(mute_posits,
                turn_to_aa(wild_seq, orf_dict[trans_id][0][0]),
                turn_to_aa(mute_seq, orf_dict[trans_id][0][0])
                )
    except:
        print "find and print kmers failure"
        raise
        return
    '''for cds_stretch in cds_list:
        (seq_start, seq_length) = cds_stretch
        try:
            wild_seq += get_seq(last_chrom, seq_start, seq_length, ref_ind)
        except:
            return
    cds_start = cds_list[0][0]
    mute_seq = make_mute_seq(wild_seq,mute_locs)
    kmer(mute_posits,
        turn_to_aa(wild_seq, orf_dict[trans_id][0][0]), 
        turn_to_aa(mute_seq, orf_dict[trans_id][0][0])
        )'''

def remove_overlaps(seq_list):
    curr_max = 0
    new_list = []
    for start_pos,stretch_length in seq_list:
        if start_pos > curr_max:
            curr_max = start_pos
        curr_end = start_pos + stretch_length
        if curr_end >= curr_max:
            if curr_max == 0:
                curr_max = start_pos
            new_list.append((curr_max, curr_end-curr_max))
            curr_max = curr_end
    new_list = list(filter(lambda x: x[1] != 0, new_list))
    return new_list

def calculate_intron_length(analysis_list):
    total_introns = 0
    for index in range(len(analysis_list)-1):
        total_introns += analysis_list[index+1][0]-(analysis_list[index][0]+analysis_list[index][1])
    return total_introns

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', type=str, required=False,
        default='-',
        help='input vcf or "-" for stdin'
    )
parser.add_argument('-x', '--bowtie-index', type=str, required=True,
        help='path to Bowtie index basename'
    )
parser.add_argument('-d', '--dicts', type=str, required=True,
        help='input path to pickled dictionaries'
    )
args = parser.parse_args()
ref_ind = bowtie_index.BowtieIndexReference(args.bowtie_index)

my_dicts = pickle.load ( open (args.dicts, "rb"))
(cds_dict, orf_dict) = (my_dicts[0], my_dicts[1])


try:
    if args.vcf == '-':
        if sys.stdin.isatty():
            raise RuntimeError('Nothing piped into this script, but input is '
                               'to be read from stdin')
        else:
            input_stream = sys.stdin
    else:
        input_stream = open(args.vcf, "r")
        last_chrom = "None" 
        line_count = 0
        seq_list = []
        #Added these vars below
        mute_seq_pos = 0
        seq_end = 0
        last_trans = 0
        total_intron_length = 0
        bounds_set = set()
        for line in input_stream:
            line_count += 1
            if not line or line[0] == '#': 
                continue
            vals = line.strip().split('\t')
            (chrom, pos, alt, info) = (vals[0], int(vals[1]), vals[4], vals[7]
                )
            tokens = info.strip().split('|')
            mute_type = tokens[1]
            if(mute_type != "missense_variant"):
                continue
            (trans_id, rel_pos) = (tokens[6], int(tokens[13]))
            if mute_type == "missense_variant":
                pos_in_codon = (rel_pos+2)%3 #ie: ATG --> 0,1,2
            try:
                if orf_dict[trans_id][0][0] == "-": 
                    pos_in_codon = 2-pos_in_codon
            except:
                #Added this line below
                trans_id = last_trans
                continue
            if last_chrom == chrom and pos <= seq_end:
                #Fixed the call with a tuple value instead of just pos
                (cds_list_left, temp, bounds_set) = get_cds(trans_id, [(pos, None)],30+pos_in_codon, 0, cds_dict, dict())
                if(len(cds_list_left)==0): continue
                adjusted_start = cds_list_left[0][0] # <-- Take out?
                (cds_list_right, temp, bounds_set) = get_cds(trans_id, [(pos, None)], 0, 32-pos_in_codon, cds_dict, dict())
                seq_list.extend(cds_list_right)
                (last_seq_start, last_seq_length) = cds_list_right.pop()
                seq_end = last_seq_start+last_seq_length
                temp_list = remove_overlaps(seq_list)
                mute_pos_in_list = bisect.bisect([pair[0] for pair in temp_list], pos) - 1
                total_intron_length = calculate_intron_length(temp_list[:mute_pos_in_list+1])
                mute_seq_pos = pos-temp_list[0][0]-total_intron_length
                #end_ind = pos+32-pos_in_codon
            else:
                if len(seq_list) != 0:
                    #(new_list, temp, bounds_set) = get_cds(last_trans, [(seq_list[0][0],None)], 0, last_pos-adjusted_start+32-last_codon_pos, cds_dict, dict())
                    new_list = remove_overlaps(seq_list)
                    find_seq_and_kmer(new_list, last_chrom, ref_ind,
                                          mute_locs, orf_dict, last_trans, mute_posits)
                (mute_locs, mute_posits, seq_list) = (dict(), [], [])
                #Fixed the call with a tuple value instead of just pos
                (cds_list, mute_locs, bounds_set) = get_cds(trans_id, [(pos, None)], 30+pos_in_codon, 32-pos_in_codon, cds_dict, mute_locs)
                if(len(cds_list)!=0):
                    seq_list.extend(cds_list)
                    (last_seq_start, last_seq_length) = cds_list.pop()
                    seq_end = last_seq_start + last_seq_length
                else:
                    chrom = "None"
                    continue
                st_ind = pos-30-pos_in_codon
                mute_seq_pos = 30+pos_in_codon
                end_ind = pos+32-pos_in_codon
            mute_locs[mute_seq_pos] = alt
            mute_posits.append((pos, line_count))
            (last_pos,last_chrom, last_trans, last_codon_pos) = (pos, chrom, trans_id, pos_in_codon)
        #Below, I changed right_side from end_ind-last_pos to seq_end
        try:
            (left_side,right_side) = (last_pos-st_ind,seq_end-last_pos)
        except NameError:
            pass
        (cds_list,mute_locs, bounds_set) = get_cds(trans_id, mute_posits, left_side, right_side, cds_dict, mute_locs)
        if(len(cds_list) != 0):
            #Below, use seq_list instead of cds_list
            new_list = remove_overlaps(seq_list)
            find_seq_and_kmer(new_list, last_chrom, ref_ind, mute_locs,
                              orf_dict, trans_id, mute_posits)

finally:
    if args.vcf != '-':
        input_stream.close()