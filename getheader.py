def getheader(vcffile):
    vcf = open(vcffile, "r")
    header = open("header.txt", "w+");
    temp = open("tempvcf.txt", "w+");
    for lines in vcf:
        if (lines[0] == '#'):
            header.write(lines)
        else:
            temp.write(lines)


getheader("testvcf.vcf")
