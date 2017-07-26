#Author: Austin Nguyen
#Date: 6-22-17
#Description: Testing parsing capabilities of pysam and also the ease of extraction/modification

import pysam
#just some test variables
'''chromosome = "NC_000017.10"
position   = 43384880
mutation   = "G"'''

paramfile = open("testinput.txt", "r")

#reads in the input file as an Alignment file
inputfile = pysam.AlignmentFile("testbam.bam", "rb") 

#fetch returns the reads with any nucleotide within the range provided
#read acts as an aligned segment
#if we used data = inputfile..., data would be an iterator which has very little documentation
#however, we can keep use "data" if we want to iterate through later.  Better modularization
#Function: Gets info from bam file and returns needed vars
#Inputs: the bam file, position of mutation, chromosome number, and the mutation base
#Outputs: Just printing and writing filtered results to bam file
#Description: Selects from bam using fetch, then matches the criteria of the correct base at the correct position
def getbaminfo(header, inputfile, position, chromosome, mutation):
    for read in inputfile.fetch(chromosome, start = position, end = position+1):
        if (read.reference_start < position and read.reference_end > position and read.query_sequence[position-read.reference_start-1] == mutation):
            print read.query_name, read.reference_start, read.reference_end, read.query_alignment_length, read.query_sequence, position-read.reference_start, read.query_sequence[position-read.reference_start-1]
            with pysam.AlignmentFile("outputbam.bam", "wb", header=header) as outf:
                a = pysam.AlignedSegment()
                a = read
                outf.write(a)

    return
#Simple Parsing of the param file, we only need 3 attributes to make the fetch from pysam

for lines in paramfile:
    parsedline = lines.strip().split()
    chromosome = parsedline[0]
    position   = int(parsedline[1])
    mutation   = parsedline[2]
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }
    getbaminfo(header, inputfile, position, chromosome, mutation)
    #print chromosome, position, mutation
    




    
