import os

def combinevcf(vcf1, vcf2, outfile="Combined.vcf"):
    vcffile = open(vcf2, "r")
    temp = open(vcf2 + ".tumortemp", "w+");
    header = open(vcf2 + ".header", "w+");
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
        else:
            header.write(lines)
    vcffile.close()
    temp.close()
    header.close()
    vcffile = open(vcf1, "r")
    temp = open(vcf2 + ".germlinetemp", "w+");
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
    vcffile.close()
    temp.close()    
    markgermline = '''awk '{print $0, "0"}' ''' + vcf2 + ".germlinetemp > " + vcf2 + ".germline"
    marktumor    = '''awk '{print $0, "1"}' ''' + vcf2 + ".tumortemp > " + vcf2 + ".tumor"
    os.system(markgermline)
    os.system(marktumor)
    command = "cat " + vcf2 + ".germline " + vcf2 + ".tumor > " + vcf2 + ".combine1"
    os.system(command)
    command2 = "sort -k1,1 -k2,2n " + vcf2 + ".combine1 > " + vcf2 + ".sorted"
    os.system(command2)
    command3 = "cat " + vcf2 + ".header " + vcf2 + ".sorted > " + vcf2 + ".combine2"
    os.system(command3)
    cut = "cut -f1,2,3,4,5,6,7,8,9,10 " + vcf2 + ".combine2 > " + outfile
    os.system(cut)
    for file in [".tumortemp", ".germlinetemp", ".combine1", ".combine2", ".sorted", ".tumor", ".germline", ".header"]:
        cleanup = "rm " + vcf2 + file
        os.system(cleanup)
