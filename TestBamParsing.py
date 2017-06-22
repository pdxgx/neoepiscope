#Author: Austin Nguyen
#Date: 6-22-17
#Description: Testing parsing capabilities of pysam and also the ease of extraction/modification

import pysam
#reads in the input file as an Alignment file
inputfile = pysam.AlignmentFile("testbam.bam", "rb") 

#fetch returns the reads with any nucleotide within the range provided
#read acts as an aligned segment
#if we used data = inputfile..., data would be an iterator which has very little documentation
#however, we can keep use "data" if we want to iterate through later.  Better modularization
for read in inputfile.fetch("NC_000017.10", start = 43384873, end = 43384986): 
    print read.reference_start, read.reference_end, read.query_alignment_length, read.query_sequence

