import pickle
import argparse
import bisect


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dump', type=str, required=True,
        help='input path to file for dict storage'
    )
parser.add_argument('-g', '--gtf', type=str, required=False,
        help='input gtf'
    )
args = parser.parse_args()

my_file = open(args.gtf, "r") # ex: open("gencode.txt").read()

cds_dict = {}
cds_orf_dict = {}
exon_dict = {}
exon_orf_dict = {}
for line in my_file:
    if not line or line[0] == '#': 
        continue
    tokens = line.strip().split('\t')
    if tokens[2] != "CDS" and tokens[2] != "exon": 
        continue
    read_info = tokens[8].split(';')
    version_in_id = read_info[1].find('.')
    if version_in_id == -1:
        transcript_id = read_info[1][16:]
    else:
        transcript_id = read_info[1][16:version_in_id]
    if transcript_id not in cds_dict and tokens[2] == "CDS":
        cds_dict[transcript_id] = [int(tokens[3]), int(tokens[4])]
        cds_orf_dict[transcript_id] = [tokens[6] + str(tokens[7])]
    elif transcript_id in cds_dict and tokens[2] == "CDS":
        insert_point = 2*bisect.bisect(cds_dict[transcript_id][0::2],
                                       int(tokens[3]))
        #I'm assuming there aren't two CDSs that start at the same point.
        cds_dict[transcript_id] = (cds_dict[transcript_id][:insert_point] 
                                    + [int(tokens[3]), int(tokens[4])] 
                                    + cds_dict[transcript_id][insert_point:])
        cds_orf_dict[transcript_id] = (cds_orf_dict[transcript_id][:insert_point//2]
                                    + [tokens[6] + str(tokens[7])]
                                    + cds_orf_dict[transcript_id][insert_point//2:])
    elif transcript_id not in exon_dict and tokens[2] == "exon":
        exon_dict[transcript_id] = [int(tokens[3]), int(tokens[4])]
        exon_orf_dict[transcript_id] = [tokens[6] + str(tokens[7])]
    elif transcript_id in exon_dict and tokens[2] == "exon":
        insert_point = 2*bisect.bisect(exon_dict[transcript_id][0::2],
                                       int(tokens[3]))
        #I'm assuming there aren't two CDSs that start at the same point.
        exon_dict[transcript_id] = (exon_dict[transcript_id][:insert_point] 
                                    + [int(tokens[3]), int(tokens[4])] 
                                    + exon_dict[transcript_id][insert_point:])
        #Exons 
        exon_orf_dict[transcript_id] = (exon_orf_dict[transcript_id][:insert_point//2]
                                    + [tokens[6] + str(tokens[7])]
                                    + exon_orf_dict[transcript_id][insert_point//2:])

    
#@TODO: Don't forget to pickle the chrom_dict also!!!!!!!
pickle_out = open(args.dump, "wb")
pickle.dump([cds_dict, cds_orf_dict, exon_dict, exon_orf_dict], pickle_out)
pickle_out.close()