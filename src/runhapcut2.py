import os

def runhapcut2(vcffile, bamfile):
        script = "HapCUT2/build/extractHAIRS --indels 1 --bam " + bamfile + " --vcf " + vcffile + " --out fragment_file"
        hap = "HapCUT2/build/HAPCUT2 --fragments fragment_file --vcf " + vcffile + " --output haplotype_output_file"

        os.system("bash -c '%s'" % script)
        os.system("bash -c '%s'" % hap)
