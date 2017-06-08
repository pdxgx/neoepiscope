#Parse vcf by creating set of tuples; form: (chrom, pos, ref, alt)
def altParseVCF(vcfStr):
	allMuts = list()
	for line in vcfStr.splitlines():
		if((len(line) < 1) or (not line[0].isdigit())):
			continue
		vals = line.strip().split('\t')
		allMuts.append((vals[0], vals[1], vals[3], vals[4]))
	return allMuts
