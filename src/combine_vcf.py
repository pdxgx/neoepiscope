import os
import getheader as gh

def combinevcf(vcf1, vcf2):
    gh.getheader(vcf2)
    markgermline = '''awk '{print $0, "0"}' ''' + vcf1 + "> germlinevcf.vcf"
    marktumor    = '''awk '{print $0, "1"}' tempvcf.txt > tempvcf.vcf'''
    os.system(markgermline)
    os.system(marktumor)
    command = "cat germlinevcf.vcf tempvcf.vcf > firstcombine.vcf"
    os.system(command)
    gh.getheader("firstcombine.vcf")
    command2 = "sort -k1,1 -k2,2n tempvcf.txt > sortedvct.txt"
    os.system(command2)
    command3 = "cat header.txt sortedvct.txt > combinedvcf.vcf"
    os.system(command3)
