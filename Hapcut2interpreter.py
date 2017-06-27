#Author: Austin Nguyen
#Date: 6-26-17
#Description: Returns phased segment of queried data.
#Input: Takes in chromosome, pos:start-end, and reference/germline sequence

#Function: returnphasing
#Inputs: string- chromsome#, int- start position, int- end position, string- reference sequence from either reference genome or germline
#Output: 1 or 2 sequences that include variants, 1 if all variants are on same chromosome, 2 if variants are spread across both chromosome
#Description: will apply mutations after searching throught the hapcut2 results
#Description: applies mutation to the relevant chromosome strand


def returnphasing(chromosome, startpos, endpos, refseq):
    chrome1 = list(refseq)
    chrome2 = list(refseq)
    chrome1count = 0 #probably unnecessary due to vcflist but good to keep counters for debug purposes
    chrome2count = 0
    chrome1vcflist = []
    chrome2vcflist = []
    errorflag = 2
    hapcutfile = open("haplotype_output_file", "r")
    for line in hapcutfile:
        stripped = line.strip().split()
        if stripped[0] != "********":
            if stripped[3]  == chromosome and int(stripped[4]) < endpos and int(stripped[4]) > startpos:
                errorflag = 0
                #print stripped[3], stripped[1], stripped[2], chromosome, stripped[4], endpos, startpos
                mutpos = int(stripped[4]) - startpos
                if int(stripped[1]) != 0:
                    chrome1[mutpos] = stripped[6]
                    chrome1count = 1
                    chrome1vcflist.append(int(stripped[0]))
                if int(stripped[2]) != 0:
                    chrome2[mutpos] = stripped[6]
                    chrome2count = 1
                    chrome2vcflist.append(int(stripped[0]))

    '''if chrome1count == 0 and chrome2count == 0:
        linecount = 0
        vcffile = open("mutect.vcf", "r")
        for line in vcffile:
            if not line or line[0] == '#': continue
            linecount += 1
            stripped = line.strip().split()
            if stripped[0] == chromosome and int(stripped[1]) < endpos and int(stripped[1]) > startpos:
                mutpos = int(stripped[1]) - startpos
                chrome1[mutpos] = stripped[4]
                chrome1[vcflist].append(linecount)
                errorflag = 1

        vcffile.close()'''
                
    hapcutfile.close()
    
    chrome1 = ''.join(chrome1)
    chrome2 = ''.join(chrome2)
    if chrome1count == 1 and chrome2count == 1:
        return (chrome1, chrome1vcflist, chrome2, chrome2vcflist, errorflag)
    if chrome1count == 1 and chrome2count == 0:
        return (chrome1, chrome1vcflist, errorflag)
    if chrome1count == 0 and chrome2count == 1:
        return (chrome2, chrome2vcflist, errorflag)
    else:
        if errorflag == 1:
            return (chrome1, chrome1vcflist, errorflag)
        else:
            return (refseq,errorflag)
            

result = returnphasing("1", 248801634, 248801640, "AAGGTT")
print result
