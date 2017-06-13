import argparse
import bowtie_index
import sys

def get_seq(start, end, chrom, ref_ind):
    splice_length = end - start 
    chr_name = "chr"+chrom #proper
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
            (chrom, pos, ref, alt) = (vals[0], int(vals[1]), vals[3], vals[4])
            if last_chrom == chrom and pos-last_pos < 33:
                #The order of that if-statement is important! Don't change it!
                end_ind = pos+32
            else:
                if last_chrom != "None":
                    seq_strand = get_seq(st_ind,end_ind,last_chrom,ref_ind)
                    #mute_strand = make_mute_strand(seq_strand,mute_locs)
                    #@TODO, now pass into makeIntoAA/ kmer function
                    #vars needed to be passed: st_ind, end_ind, last_chrom,
                    #seq_strand, mute_strand
                mute_locs = dict()
                (st_ind, end_ind) = (pos-32, pos+32)
            mute_locs[(pos-st_ind)] = alt
            (last_pos,last_chrom) = (pos, chrom)
        seq_strand = get_seq(st_ind, end_ind, last_chrom, ref_ind)
        #mute_strand = make_mute_strand(seq_strand,mute_locs)
        #@TODO, now pass into makeIntoAA/ kmer function
        #vars needed to be passed: st_ind, end_ind, last_chrom,
        #seq_strand, mute_strand
    #@TODO Repeated code above; need to clean/ make helper function
finally:
    if args.vcf != '-':
        input_stream.close()
