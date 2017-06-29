startpos = 100
endpos = 114

refseq    = "AGTGCGTGGTGTGTG"

queryseq= "ACG"
queryseq2= "AAGTGGTTGGTGTGTGG"
queryseq3= "AAGTGGT"
queryseq4= "GCGG"
querystart = 103
queryend = 105
querystart2 = 99
queryend2 = 115
querystart3 = 99
queryend3 = 105
querystart4 = 112
queryend4 = 115
strand1 = list(refseq)
passed1 = [0] * len(refseq)
mutlist1 = []
mutpos1 = []
strand2 = list(refseq)
passed2 = [0] * len(refseq)
mutlist2 = []
mutpos2 = []

def testfunc(queryseq, querystart, queryend):
    if querystart > startpos and queryend < endpos: #case 2
        counter = querystart-startpos
        applyseq = 1
        for base in queryseq:
            if passed1[counter] == 0:
                        counter += 1
            else:
                if strand1[counter] == base:
                    counter += 1
                else:
                    applyseq = 2
                    break
        counter = querystart-startpos
        if applyseq == 1:
            for base in queryseq:
                if strand1[counter] != base:
                    mutlist1.append(base)
                    mutpos1.append(counter+startpos)
                strand1[counter] = base
                passed1[counter] = 1
                counter += 1
                
        else:
            for base in queryseq:
                if strand2[counter] != base:
                   mutlist2.append(base)
                   mutpos2.append(counter+startpos)           
                strand2[counter] = base
                passed2[counter] = 1
                counter += 1

    if querystart < startpos and queryend > endpos: #case 1
        counter = 0
        applyseq = 1
        for x in range(0, abs(querystart-startpos)):
            queryseq = queryseq[1:]

        for base in queryseq:
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
            for base in queryseq:
                if strand1[counter] != base:
                    mutlist1.append(base)
                    mutpos1.append(counter+startpos)
                strand1[counter] = base
                passed1[counter] = 1
                counter += 1
                
        else:
            for base in queryseq:
                if strand2[counter] != base:
                   mutlist2.append(base)
                   mutpos2.append(counter+startpos)           
                strand2[counter] = base
                passed2[counter] = 1
                counter += 1

    if querystart < startpos and queryend < endpos: #case 3
        counter = 0
        applyseq = 1
        for x in range(0, abs(querystart-startpos)):
            queryseq = queryseq[1:]
        for base in queryseq:
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
            for base in queryseq:
                if strand1[counter] != base:
                    mutlist1.append(base)
                    mutpos1.append(counter+startpos)
                strand1[counter] = base
                passed1[counter] = 1
                counter += 1
                
        else:
            for base in queryseq:
                if strand2[counter] != base:
                   mutlist2.append(base)
                   mutpos2.append(counter+startpos)           
                strand2[counter] = base
                passed2[counter] = 1
                counter += 1

    if querystart > startpos and queryend > endpos: #case4
        counter = querystart-startpos
        applyseq = 1
        for x in range(0, (queryend-endpos)):
            queryseq = queryseq[:-1]
        for base in queryseq:
            if passed1[counter] == 0:
                        counter += 1
            else:
                if strand1[counter] == base:
                    counter += 1
                else:
                    applyseq = 2
                    break
        counter = querystart-startpos
        if applyseq == 1:
            for base in queryseq:
                if strand1[counter] != base:
                    mutlist1.append(base)
                    mutpos1.append(counter+startpos)
                strand1[counter] = base
                passed1[counter] = 1
                counter += 1
                
        else:
            for base in queryseq:
                if strand2[counter] != base:
                   mutlist2.append(base)
                   mutpos2.append(counter+startpos)           
                strand2[counter] = base
                passed2[counter] = 1
                counter += 1

                
    print ''.join(strand1)
    print mutlist1
    print mutpos1
    print passed1
    print ''.join(strand2)
    print mutlist2
    print mutpos2
    print passed2

print "Seq 1"
testfunc(queryseq, querystart, queryend)
#print "Seq 2"
#testfunc(queryseq2, querystart2, queryend2)
#print "Seq 3"
#testfunc(queryseq3, querystart3, queryend3)
print "Seq 4"
testfunc(queryseq4, querystart4, queryend4)
