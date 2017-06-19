import pickle
import argparse
import bisect


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dump', type=str, required=True,
        help='input path to file for dict storage'
    )

args = parser.parse_args()

my_file = ""# ex: open("gencode.txt").read()

exon_dict = {}
chrom_dict = {}
for line in my_file.splitlines():
    if not line or line[0] == '#': continue
    tokens = line.strip().split('\t')
    if tokens[2] != "exon": continue
    read_info = tokens[8].split(';')
    version_in_id = read_info[1].find('.')
    if version_in_id == -1:
        transcript_id = read_info[1][16:]
    else:
        transcript_id = read_info[1][16:version_in_id]
    if transcript_id not in exon_dict:
        exon_dict[transcript_id] = [int(tokens[3]), int(tokens[4])]
    else:
        insert_point = 2*bisect.bisect(exon_dict[transcript_id][0::2],
                                       int(tokens[3]))
        #I'm assuming there aren't two exons that start at the same point.
        exon_dict[transcript_id] = (exon_dict[transcript_id][:insert_point] 
                                    + [int(tokens[3]), int(tokens[4])] 
                                    + exon_dict[transcript_id][insert_point:])
    if transcript_id not in chrom_dict:
        chrom_dict[transcript_id] = tokens[0]
    
#@TODO: Don't forget to pickle the chrom_dict also!!!!!!!
pickle_out = open(args.dump, "wb")
pickle.dump(exon_dict, pickle_out)
pickle_out.close()