#Summary: Tumor filter will take results from previous blackbox in pipeline and eliminate reads that do not contain those mutated strings
#Input: Tumor bam/sam file, mutated nucleotide k-mers
#Process: Filter out reads from bam/sam that do not have the mutated nucleotide k-mers in alignment
#Output: Only tumor reads that have mutations in list given

