import os

def combinevcf(vcf1, vcf2, outfile="Combined.vcf"):
    vcffile = open(vcf2, "r")
    temp = open(vcf2 + ".temp", "w+");
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
    vcffile.close()
    temp.close()
    markgermline = '''awk '{print $0, "0"}' ''' + vcf1 + "> " + vcf2 + ".germline"
    marktumor    = '''awk '{print $0, "1"}' ''' + vcf2 + ".temp > " + vcf2 + ".temp2"
    os.system(markgermline)
    os.system(marktumor)
    command = "cat " + vcf2 + ".germline " + vcf2 + ".temp2 > " + vcf2 + ".combine1"
    os.system(command)
    vcffile = open(vcf2 + ".combine1", "r")
    header = open(vcf2 + ".header", "w+");
    for lines in vcffile:
        if (lines[0] == '#'):
            header.write(lines)
    header.close()
    command2 = "sort -k1,1 -k2,2n " + vcf2 + ".temp2 > " + vcf2 + ".sorted"
    os.system(command2)
    command3 = "cat " + vcf2 + ".header " + vcf2 + ".sorted > " + vcf2 + ".combine2"
    os.system(command3)
    cut = "cut -f1,2,3,4,5,6,7,8,9,10 " + vcf2 + ".combine2 > " + outfile
    os.system(cut)
#    for file in [".temp", ".temp2", ".combine2", ".combine1", ".sorted", ".germline", ".header"]:
 #       cleanup = "rm " + vcf2 + file
  #      os.system(cleanup)
