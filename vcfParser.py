import argparse
import bowtie_index
import sys

def get_seq(start, end, chrom, ref_ind):
    splice_length = end - start 
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
                    seq_strand = get_seq(st_ind,end_ind,last_chrom,ref_ind)
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
