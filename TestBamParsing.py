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

def getbaminfo(inputfile, position, chromosome, mutation):
    for read in inputfile.fetch(chromosome, start = position, end = position+1):
        if (read.reference_start < position and read.reference_end > position and read.query_sequence[position-read.reference_start-1] == mutation):
            print read.reference_start, read.reference_end, read.query_alignment_length, read.query_sequence, position-read.reference_start, read.query_sequence[position-read.reference_start-1]
            #print read

    return
#Simple Parsing of the param file, we only need 3 attributes to make the fetch from pysam

for lines in paramfile:
    parsedline = lines.strip().split()
    chromosome = parsedline[0]
    position   = int(parsedline[1])
    mutation   = parsedline[2]
    getbaminfo(inputfile, position, chromosome, mutation)
    #print chromosome, position, mutation
    

#fetch returns the reads with any nucleotide within the range provided
#read acts as an aligned segment
#if we used data = inputfile..., data would be an iterator which has very little documentation
#however, we can keep use "data" if we want to iterate through later.  Better modularization



    
