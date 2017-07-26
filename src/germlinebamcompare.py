#Author: Austin Nguyen
#Date: 6-28-17
#Description: Takes reference query and compares it to the germline bam
#Returns: one or more sequences of germline bam

#Function: compares ref sequence to germline bam,
#returns at most 2 possible sequences if mutations have differences
#Works by passing through strand1 for every read. 
#if strand has already been passed once and still differences, applies changes to strand 2
#inputs: chromosome, startpos, endpos, reference sequence
#outputs: 1 or 2 strands depending on where the read mutations are
#
def comparegermline(chromosome, startpos, endpos, refseq):
    germlinefile = pysam.AlignmentFile("germlinebam.bam", "rb")
    strand1 = list(refseq)
    passed1 = [0] * len(refseq)
    mutlist1 = []
    mutpos1 = []
    strand2 = list(refseq)
    passed2 = [0] * len(refseq)
    mutlist2 = []
    mutpos2 = []
    for read in germlinefile.fetch(chromosome, start = startpos, end = endpos):
        if read.reference_start < startpos and read.reference_end > endpos: #case 1
            counter = 0
            applyseq = 1
            for x in range(0, abs(read.reference_start - startpos)):
                read.query_sequence[1:]
            for x in range(0, (read.reference_end - endpos)):
                read.query_sequence[:-1]
            for base in read.query_sequence:
                if passed1[counter] == 0:
                    counter += 1
                else:
                    if strand1[counter] == base:
                        counter += 1
                    else:
                        applyseq = 2
                        break

            counter = 0
            if applyseq == 1:
                for base in read.query_sequence:
                    if strand1[counter] != base:
                        mutlist1.append(base)
                        mutpos1.append(counter+startpos)
                    strand1[counter] = base
                    passed2[counter] = 1
                    counter += 1

            else:
                for base in read.query_sequence:
                    if strand2[counter] != base:
                        mutlist2.append(base)
                        mutpos2.append(counter+startpos)
                    strand2[counter] = base
                    passed2[counter] = 1
                    counter += 1

        if read.reference_start > startpos and read.reference_end < endpos: #case 2
            counter = read.reference_start-startpos
            applyseq = 1
            for base in read.query_sequence:
                if passed1[counter] == 0:
                    counter += 1
                else:
                    if strand1[counter] == base:
                        counter += 1
                    else:
                        applyseq = 2
                        break

            counter = read.reference_start-startpos
            if applyseq == 1:
                for base in read.query_sequence:
                    if strand1[counter] != base:
                        mutlist1.append(base)
                        mutpos1.append(counter+startpos)
                    strand1[counter] = base
                    passed2[counter] = 1
                    counter += 1
                    
            else:
                for base in read.query_sequence:
                    if strand2[counter] != base:
                        mutlist2.append(base)
                        mutpos2.append(counter+startpos)
                    strand2[counter] = base
                    passed2[counter] = 1
                    counter += 1
                    

        if read.reference_start < startpos and read.reference_end < endpos: #case 3
            counter = 0
            applyseq = 1
            for x in range(0, abs(read.reference_start - startpos)):
                read.query_sequence[1:]
            for base in read.query_sequence:
                if passed1[counter] == 0:
                    counter += 1
                else:
                    if strand1[counter] == base:
                        counter += 1
                    else:
                        applyseq = 2
                        break

            counter = 0
            if applyseq == 1:
                for base in read.query_sequence:
                    if strand1[counter] != base:
                        mutlist1.append(base)
                        mutpos1.append(counter+startpos)
                    strand1[counter] = base
                    passed2[counter] = 1
                    counter += 1
                    
            else:
                for base in read.query_sequence:
                    if strand2[counter] != base:
                        mutlist2.append(base)
                        mutpos2.append(counter+startpos)
                    strand2[counter] = base
                    passed2[counter] = 1
                    counter += 1
        if read.reference_start > startpos and read.reference_end > endpos: #case 4
            counter = read.reference_start-startpos
            applyseq = 1
            for x in range(0, (read.reference_end - endpos)):
                read.query_sequence[:-1]
            for base in read.query_sequence:
                if passed1[counter] == 0:
                    counter += 1
                else:
                    if strand1[counter] == base:
                        counter += 1
                    else:
                        applyseq = 2
                        break

            counter = read.reference_start-startpos
            if applyseq == 1:
                for base in read.query_sequence:
                    if strand1[counter] != base:
                        mutlist1.append(base)
                        mutpos1.append(counter+startpos)
                    strand1[counter] = base
                    passed2[counter] = 1
                    counter += 1
                    
            else:
                for base in read.query_sequence:
                    if strand2[counter] != base:
                        mutlist2.append(base)
                        mutpos2.append(counter+startpos)
                    strand2[counter] = base
                    passed2[counter] = 1
                    counter += 1
    
        
    
