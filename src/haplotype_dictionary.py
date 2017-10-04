import collections
from sortedcontainers import SortedDict
import re

def getfreqlabel(vcffile):
    vcf = open(vcffile, "r")
    for read in vcf:
        flag = 1
        if 'allele frequency' in read.lower():
                freqlabel = re.sub(r'.*ID=([A-Z]+)[,].*', r'\1', read)
                '''if 'info' in read.lower():
                    startchar = 11
                else:
                    startchar = 13
                freqlabel = ""
                while read[startchar] != ',':
                    freqlabel = freqlabel + read[startchar] 
                    startchar += 1'''
                #print freqlabel
                flag = 0
                break
    vcf.close()
    if flag == 1:
        print "Raise runtime error"
    return freqlabel

def create_haplotype_dictionary(haplooutfile, freqpos):
    '''
        Returns dictionary of haplotype information
        haplooutfile: name of haplotype output file
        freqpos: location in genotype field where allele frequency is (0 based)

        outputs: dictionary with key of (chrome, pos) and value of (allele on chrome A, allele on chrome B, position, reference, variant, allele frequency, phasing score)
    '''
    hapfile = open(haplooutfile, "r")
    haplotype_dict = collections.OrderedDict()
    for line in hapfile:
        hapline = line.strip().split()
        if hapline[0] != "********" and hapline[0] != "BLOCK:":
            infoline = hapline[7].strip().split(':')
            haplotype_dict[(int(hapline[3]), int(hapline[4]))] = (int(hapline[1]), int(hapline[2]), hapline[5], hapline[6], float(infoline[freqpos].rstrip('%')), float(hapline[10]))
    return SortedDict(haplotype_dict)


blah = create_haplotype_dictionary("haplotype_temp", 5)
for k, v in blah.items():
    print k, v
