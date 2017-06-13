import argparse
import bowtie_index
import sys

def get_seq(start, end, chrom, ref_ind):
    splice_length = end - start
    return(ref_ind.get_stretch(chrom, start, splice_length))

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
        #@TODO Try-except in case invalid file type
        last_pos = -1
        info_list = []
        mute_locs = set()
        for line in input_stream:
            line = line.strip()
            if not line or line[0] == '#': continue
            vals = line.strip().split('\t')
            (chrom, pos, ref, alt) = (vals[0], int(vals[1]), vals[3], vals[4])
            if last_pos == -1:
                (st_ind, end_ind, last_pos) = (pos-32, pos+32, pos)
                last_chrom = chrom
                mute_locs.add((pos - st_ind, alt)) #Establish base case
            else:
                if pos-last_pos < 33 or last_chrom != chrom:
                    end_ind = pos+32
                else:
                    info_list.append((st_ind,end_ind,mute_locs))
                    mute_locs = set()
                    (st_ind, end_ind) = (pos-32, pos+32)
                mute_locs.add((pos-st_ind,alt))
                last_pos = pos
            info_list.append((st_ind,end_ind,mute_locs))
            seq_splice = get_seq(st_ind,end_ind,last_chrom,ref_ind)
            print(seq_splice)
finally:
    if args.vcf != '-':
        input_stream.close()
