import pickle
import re

# bisect_left function
#interval tree interval 1.0.0 library

#define coordinate input --> transcript label
# transcript label input --> cds output

def gtf_to_cds(gtf_file, pickle_dict = ""):
        gtf_data = open(gtf_file, "r")
        #cds_dict[transcript_id]= [[start, stop, 1(if CDS)/0(if stop), reading frame, +/-], [start, stop, 1(if CDS)/0(if stop), reading frame, +/-], ...]
        cds_dict = {}
        relevant_identifiers = set(["CDS", "start_codon", "stop_codon"])
        for line in gtf_data:
                if not line or line[0] == '#':
                        continue
                tokens = line.strip().split('\t')
                if (tokens[2] != 'CDS' and tokens[2] != 'stop_codon'):
                        continue    
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

        # sort cds_dict coordinates (left -> right) for each transcript                                
        for transcript_id in cds_dict:
                cds_dict[transcript_id].sort(key=lambda x: x[0])
        gtf_data.close()
        
        if pickle_dict != "":
                pickle_out = open(pickle_dict, "wb")
                pickle.dump(cds_dict, pickle_out)
                pickle_out.close()    
                
        return cds_dict



#    with open(args.dicts, 'rb') as dict_stream:
 #       (cds_dict) = pickle.load(dict_stream)
