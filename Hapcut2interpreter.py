#Author: Austin Nguyen
#Date: 6-26-17
#Description: Returns phased segment of queried data.
#Input: Takes in chromosome, pos:start-end, and reference/germline sequence

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
    hapscore1 = [0] * len(refseq)
    hapscore2 = [0] * len(refseq)
    chromepassed2 = [0] * len(refseq)
    chrome1count = 0 #probably unnecessary due to vcflist but good to keep counters for debug purposes
    chrome2count = 0
    chrome1vcflist = []
    chrome2vcflist = []
    errorflag = 2
    hapcutfile = open("haplotype_output_file", "r")
    for line in hapcutfile:
        stripped = line.strip().split()
        if stripped[0] != "********":
            if stripped[3]  == chromosome and int(stripped[4]) <= endpos and int(stripped[4]) >= startpos:
                errorflag = 0
                #print stripped[3], stripped[1], stripped[2], chromosome, stripped[4], endpos, startpos
                mutpos = int(stripped[4]) - startpos
                if int(stripped[1]) != 0:
                    chrome1[mutpos] = stripped[6]
                    chrome1count = 1
                    chrome1vcflist.append(int(stripped[0]))
                    chromepassed1[mutpos] = 1
                    hapscore1[mutpos] = stripped[10]
                if int(stripped[2]) != 0:
                    chrome2[mutpos] = stripped[6]
                    chrome2count = 1
                    chrome2vcflist.append(int(stripped[0]))
                    chromepassed2[mutpos] = 1
                    hapscore2[mutpos] = stripped[10]

    linecount = 0
    vcffile = open(vcfname, "r")
    for line in vcffile:
        #print line
        if not line or line[0] == '#': continue
        linecount += 1
        stripped = line.strip().split()
        if stripped[0] == chromosome and int(stripped[1]) <= endpos and int(stripped[1]) >= startpos:
            mutpos = int(stripped[1]) - startpos
            if chromepassed1[mutpos] != 1:
                chrome1[mutpos] = stripped[4]
                chrome1vcflist.append(linecount)
                errorflag = 1
                chromepassed1[mutpos] = 2
            if chromepassed2[mutpos] != 1:
                chrome2[mutpos] = stripped[4]
                chrome2vcflist.append(linecount)
                errorflag = 1
                chromepassed2[mutpos] = 2

    vcffile.close()
                   
    hapcutfile.close()
    chrome1 = ''.join(chrome1)
    chrome2 = ''.join(chrome2)
    if chrome1count == 1 and chrome2count == 1:
        return (chrome1, chrome1vcflist, chrome2, chrome2vcflist, errorflag, chromepassed1, chromepassed2, hapscore1, hapscore2)
    if chrome1count == 1 and chrome2count == 0:
        return (chrome1, chrome1vcflist, errorflag, chromepassed1, hapscore1)
    if chrome1count == 0 and chrome2count == 1:
        return (chrome2, chrome2vcflist, errorflag, chromepassed2, hapscore2)
    else:
        return (refseq,errorflag)
            

