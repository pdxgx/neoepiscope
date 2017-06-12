import argparse, bowtie_index

parser = argparse.ArgumentParser()
parser.add_argument('inp')
args = parser.parse_args()

def getSeq(start, end, chrom):
    ref_ind = bowtie_index.BowtieIndexReference('/mnt/c/Users/mihir.DESKTOP-QKPMQOM/Documents/OHSU/hg19')
    splice_length = end-start
    return(reference_index.get_stretch(chrom, start, splice_length))


with open(args.inp) as input_stream:
    #@TODO Try-except in case invalid file type
    lastPos = -1
    info_list = list()
    mute_locs = set()
    for line in input_stream:
        if((line[0]=="#") or (len(line)==0)): continue
        vals = line.strip().split('\t')
        (chrom, pos, ref, alt) = (vals[0], int(vals[1]), vals[3], vals[4])
        if(lastPos == -1):
            (st_ind, end_ind, last_pos) = (pos-32, pos+32, pos)
            mute_locs.add((pos-st_ind,alt))
        else:
            if(pos-last_pos<33): end_ind = pos+32
            else:
                info_list.append((st_ind,end_ind,mute_locs))
                mute_locs = set()
                (st_ind, end_ind) = (pos-32, pos+32)
            mute_locs.add((pos-st_ind,alt))
            lastPos = pos
        info_list.append((st_ind,end_ind,mute_locs))
        seq_splice = getSeq(st_ind,end_ind,vals[0])
        print(seq_splice)