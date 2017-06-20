import bisect
import argparse
import bowtie_index
import sys
#import pickle to unpickle ordered_exon_dict

def get_exons(transcript_id, mutation_pos, strand_length_left, 
              strand_length_right):
    ordered_exon_dict = {}
    total_strand_length = strand_length_right + strand_length_left
    original_length_left = strand_length_left
    exon_list = ordered_exon_dict[transcript_id]
    middle_exon_index = 2*bisect.bisect(exon_list[::2], mutation_pos)-2
    print(exon_list[middle_exon_index])
    #If the middle_exon_index is past the last boundary, move it to the last.
    if middle_exon_index > len(exon_list)-1:
        middle_exon_index -= 2
    nucleotide_index_list = []
    curr_left_index = middle_exon_index
    curr_right_index = middle_exon_index+1 #Exon boundary end indexes
    curr_pos_left = mutation_pos
    curr_pos_right = mutation_pos #Actual number in chromosome
    #If the mutation is not on in exon bounds, return [].
    if mutation_pos > exon_list[curr_right_index]:
        return nucleotide_index_list
    count = 0
    while(len(nucleotide_index_list) == 0 or 
          sum([index[1] for index in nucleotide_index_list]) 
          < (original_length_left)):
        if curr_pos_left-exon_list[curr_left_index] >= strand_length_left:
            nucleotide_index_list.append((curr_pos_left-strand_length_left,
                                          strand_length_left))
            strand_length_left = 0
        else:
            nucleotide_index_list.append((exon_list[curr_left_index],
                                    curr_pos_left-exon_list[curr_left_index]))
            strand_length_left -= curr_pos_left-exon_list[curr_left_index]
            curr_pos_left = exon_list[curr_left_index-1]
            curr_left_index -= 2
            if curr_left_index < 0:
                print("Exceeded all possible exon boundaries!")
                #Changed total_strand_length for comparison in next while loop.
                total_strand_length = (original_length_left
                                      - strand_length_left
                                      + strand_length_right)
                break
    while(len(nucleotide_index_list) == 0 or 
              sum([index[1] for index in nucleotide_index_list]) 
              < (total_strand_length)):
        if exon_list[curr_right_index] >= curr_pos_right + strand_length_right:
            nucleotide_index_list.append((curr_pos_right, strand_length_right))
            strand_length_right = 0
        else:
            try:
                nucleotide_index_list.append((curr_pos_right+1,
                                              exon_list[curr_right_index]
                                              - curr_pos_right))
                strand_length_right -= exon_list[curr_right_index]-curr_pos_right
                curr_pos_right = exon_list[curr_right_index+1]
                curr_right_index += 2
            except IndexError:
                print("Exceeded all possible exon boundaries!")
                break
    return nucleotide_index_list

def get_seq(chrom, start, splice_length, ref_ind):
    chr_name = "chr"+chrom #proper
    print(start, end, chr_name)
    try:
        strand = ref_ind.get_stretch(chr_name, start, splice_length)
        return strand
    except:
        return "No"

def make_mute_strand(orig_strand, mute_locs):
    mute_strand = ""
    for ind in range(len(orig_strand)):
        if ind in mute_locs:
            mute_strand += mute_locs[ind]
        else:
            mute_strand += orig_strand[ind]
    return mute_strand

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', type=str, required=False,
        default='-',
        help='input vcf or "-" for stdin'
    )
parser.add_argument('-x', '--bowtie-index', type=str, required=True,
        help='path to Bowtie index basename'
    )
args = parser.parse_args()
ref_ind = bowtie_index.BowtieIndexReference(args.bowtie_index)

try:
    if args.vcf == '-':
        if sys.stdin.isatty():
            raise RuntimeError('Nothing piped into this script, but input is '
                               'to be read from stdin')
        else:
            input_stream = sys.stdin
    else:
        input_stream = open(args.vcf)
        last_chrom = "None" #Will this work?
        for line in input_stream:
            if not line or line[0] == '#': continue
            vals = line.strip().split('\t')
            (chrom, pos, alt, info) = (vals[0], int(vals[1]), vals[4], vals[7]
                )
            tokens = info.strip().split('|')
            mute_type = tokens[1]
            if(mute_type != "missense_variant"): continue
            (trans_id, rel_pos) = (tokens[6], int(tokens[13]))
            pos_in_codon = (rel_pos+2)%3 #ATG --> 0,1,2
            if last_chrom == chrom and pos-last_pos < (32-pos_in_codon):
                #Does it matter if mutations on same transcript?
                #The order of that if-statement is important! Don't change it!
                end_ind = pos+32-pos_in_codon
            else:
                if last_chrom != "None":
                    (left_side,right_side) = (last_pos-st_ind,end_ind-last_pos)
                    exon_list = get_exons(trans_id, last_pos, left_side, right_side)
                    seq_strand = ""
                    for exon_stretch in exon_list:
                        (seq_start, strand_length) = exon_stretch
                        seq_strand += get_seq(last_chrom, seq_start, strand_length, ref_ind)
                    mute_strand = make_mute_strand(seq_strand,mute_locs)
                    #@TODO, now pass into makeIntoAA/ kmer function
                    #vars needed to be passed: st_ind, end_ind, last_chrom,
                    #seq_strand, mute_strand
                mute_locs = dict()
                st_ind = pos-30-pos_in_codon
                end_ind = pos+33-pos_in_codon
            mute_locs[(pos-st_ind)] = alt
            (last_pos,last_chrom) = (pos, chrom)
        seq_strand = get_seq(st_ind, end_ind, last_chrom, ref_ind)
        mute_strand = make_mute_strand(seq_strand,mute_locs)
        #@TODO, now pass into makeIntoAA/ kmer function
        #vars needed to be passed: st_ind, end_ind, last_chrom,
        #seq_strand, mute_strand
    #@TODO Repeated code above; need to clean/ make helper function
finally:
    if args.vcf != '-':
        input_stream.close()