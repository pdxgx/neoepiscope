import pickle
import argparse
import re


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dump', type=str, required=True,
        help='input path to file for dict storage'
    )
parser.add_argument('-g', '--gtf', type=str, required=False,
        help='input gtf'
    )
args = parser.parse_args()

my_file = open(args.gtf, "r") # ex: open("gencode.txt").read()

def gtf_to_cds(gtf_file, write_pickle_dict = False):
        gtf_data = open(gtf_file, "r")
        cds_dict = {}
        relevant_identifiers = set(["CDS", "start_codon", "stop_codon"])
        for line in gtf_data:
                if not line or line[0] == '#':
                        continue
                tokens = line.strip().split('\t')
                if (tokens[2] != 'CDS' and tokens[2] != 'stop_codon'):
                        continue                
        return cds_dict

#cds_dict[transcript_id]= [[start, stop, 1(if CDS)/0(if stop), reading frame, +/-], [start, stop, 1(if CDS)/0(if stop), reading frame, +/-], ...]
cds_dict = {}
for line in my_file:
    if not line or line[0] == '#': 
        continue
    tokens = line.strip().split('\t')
    if tokens[2] not in ["CDS", "start_codon", "stop_codon"]: 
        continue
    ## not sure if we want to be using the start_codon at all or if this info is entirely redundant (?)
    transcript_id = re.sub(r'.*transcript_id \"([A-Z0-9._]+)\"[;].*', r'\1', tokens[8])
    if transcript_id not in cds_dict:
        if tokens[2] == "CDS":
                cds_dict[transcript_id] = [[int(tokens[3]), int(tokens[4]), 1, tokens[6], int(tokens[7])]]
        elif tokens[2] == "stop_codon":
                cds_dict[transcript_id] = [[int(tokens[3]), int(tokens[4]), 0, tokens[6], int(tokens[7])]]
    else:
        if tokens[2] == "CDS":
                cds_dict[transcript_id] = cds_dict[transcript_id].append([int(tokens[3]), int(tokens[4]), 1, tokens[6], int(tokens[7])])
        elif tokens[2] == "stop_codon":
                cds_dict[transcript_id] = cds_dict[transcript_id].append([int(tokens[3]), int(tokens[4]), 0, tokens[6], int(tokens[7])])

for transcript_id in cds_dict:
        cds_dict[transcript_id].sort(key=lambda x: x[0])
                
    
pickle_out = open(args.dump, "wb")
pickle.dump(cds_dict, pickle_out)
pickle_out.close()



#    with open(args.dicts, 'rb') as dict_stream:
 #       (cds_dict) = pickle.load(dict_stream)
