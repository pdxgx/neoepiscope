import os
import getheader as gh

def combinevcf(vcf1, vcf2):
    gh.getheader(vcf2)
    command = "cat " + vcf1 + " tempvcf.txt > firstcombine.vcf"
    os.system(command)
    gh.getheader("firstcombine.vcf")
    command2 = "sort -k1,1 -k2,2n tempvcf.txt > sortedvct.txt"
    os.system(command2)
    command3 = "cat header.txt sortedvct.txt > combinedvcf.vcf"
    os.system(command3)
