#Author: Austin Nguyen

#Function: gets header and splits vcffile into 2 txt files (mutation lines and header)
#textfiles can be read back once you cat tempvcfs together or do this for one file and cat the other with it

def getheader(vcffile):
    vcf = open(vcffile, "r")
    header = open("header.txt", "w+");
    temp = open("tempvcf.txt", "w+");
    for lines in vcf:
        if (lines[0] == '#'):
            header.write(lines)
        else:
            temp.write(lines)

