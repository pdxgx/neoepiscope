import os

def runhapcut2(vcffile, bamfile, outfile="Haplotype_output_file"):
        script = "HapCUT2/build/extractHAIRS --indels 1 --bam " + bamfile + " --vcf " + vcffile + " --out " + outfile + ".frag"
        hap = "HapCUT2/build/HAPCUT2 --fragments " + outfile + ".frag --vcf " + vcffile + " --output " + outfile
        cleanup = "rm " + outfile + ".frag"

        os.system("bash -c '%s'" % script)
        os.system("bash -c '%s'" % hap)
        os.system("bash -c '%s'" % cleanup)
