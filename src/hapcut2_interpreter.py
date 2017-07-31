#Author: Austin Nguyen
#Date: 6-26-17
#Description: Returns phased segment of queried data.
#Input: Takes in chromosome, pos:start-end, and reference/germline sequence

        


#Function: get freq
#Inputs: vcf file
#Output: the position of the allele frequency from the VCF 10th column (which pos)
#Description: will read to header file and pattern match on the line

#pattern match beginning of header file VCF
#Note: check if all vcf have the the header file format
#if not, we're in trouble....
def getfreqlabel(vcffile):
    vcf = open(vcffile, "r")
    for read in vcf:
        flag = 1
        if 'allele frequency' in read.lower():
                startchar = 11
                freqlabel = ""
                while read[startchar] != ',':
                    freqlabel = freqlabel + read[startchar] 
                    startchar += 1
                #print freqlabel
                flag = 0
                break

    if flag == 1:
        print "Raise runtime error"
    return freqlabel


#Function: returnphasing
#Inputs: string- chromsome#, int- start position, int- end position, string- reference sequence from either reference genome or germline
#Output: 1 or 2 sequences that include variants, 1 if all variants are on same chromosome, 2 if variants are spread across both chromosome
#Description: will apply mutations after searching throught the hapcut2 results
#Description: will apply mutations of both hapcut2 results and vcf file (nonphased mutations)
#Description: applies mutation to the relevant chromosome strand

#Return flags: 0- mutation from Hapcut2 output, 1- mutation from vcf, 2- no mutation

def returnphasing(chromosome, startpos, endpos, refseq, vcfname):
    chrome1 = list(refseq)
    chrome2 = list(refseq)
    chromepassed1 = [0] * len(refseq)
    chromepassed2 = [0] * len(refseq)
    hapscore1 = [0] * len(refseq)
    hapscore2 = [0] * len(refseq)
    allelefreq = [0] * len(refseq)
    allelefreq = [0] * len(refseq)
    chrome1count = 0 #probably unnecessary due to vcflist but good to keep counters for debug purposes
    chrome2count = 0
    chrome1vcflist = []
    chrome2vcflist = []
    chromeflag1 = 0
    chromeflag2 = 0
    hapflag = 0
    errorflag = 0
    hapcutfile = open("haplotype_output_file", "r")
    for line in hapcutfile:
        stripped = line.strip().split()
        if stripped[0] != "********":
            if stripped[3]  == chromosome and int(stripped[4]) <= endpos and int(stripped[4]) >= startpos:
                #print stripped[3], stripped[1], stripped[2], chromosome, stripped[4], endpos, startpos
                mutpos = int(stripped[4]) - startpos
                if (stripped[1] == '-' or stripped[2] == '-'):
                    continue
                if int(stripped[1]) != 0:
                    hapflag = 1
                    chromeflag1 = 1
                    chrome1[mutpos] = stripped[6]
                    if len(stripped[5]) > 1:
                        for x in range(0,len(stripped[5])-1):
                            chrome1[mutpos+x] = ""
                    chrome1count = 1
                    chrome1vcflist.append(int(stripped[0]))
                    chromepassed1[mutpos] = 1
                    hapscore1[mutpos] = float(stripped[10])
                if int(stripped[2]) != 0:
                    hapflag = 1
                    chromeflag2 = 1
                    chrome2[mutpos] = stripped[6]
                    if len(stripped[5]) > 1:
                        for x in range(0,len(stripped[5])-1):
                            chrome1[mutpos+x+1] = ""
                    chrome2count = 1
                    chrome2vcflist.append(int(stripped[0]))
                    chromepassed2[mutpos] = 1
                    hapscore2[mutpos] = float(stripped[10])
                   
    hapcutfile.close()
    #chrome1 = ''.join(chrome1)
    #chrome2 = ''.join(chrome2)
    seqlist = []
    seqset = []
    seqlist.append(chrome1)
    seqlist.append(chrome2)
    linecount = 0
    vcffile = open(vcfname, "r")
    freqlabel = getfreqlabel(vcfname)
    for line in vcffile:
        if not line or line[0] == '#': continue
        linecount += 1
        infostring = []
        stripped = line.strip().split()
        info = stripped[8].split(':')
        allelefreqline = stripped[9].split(':')
        if stripped[0] == chromosome and int(stripped[1]) <= endpos and int(stripped[1]) >= startpos:
            mutpos = int(stripped[1]) - startpos
            if chromepassed1[mutpos] != 1 and chromepassed2[mutpos] != 1:
                length = len(seqlist)
                for x in range(0,length):
                    newseq = list(seqlist[x])
                    newseq[mutpos] = stripped[4]
                    if len(stripped[3]) > 1:
                        for x in range(0, len(stripped[3])-1):
                            newseq[mutpos+x+1] = ""
                    #newseq = ''.join(newseq)
                    seqlist.append(newseq)
    for seq in seqlist:
        seqset.append("".join(seq))
    seqset = set(seqset)
    seqback = list(seqset)
    #seqset = set(seqlist)
    #seqback = list(seqset)
    return seqback

            

