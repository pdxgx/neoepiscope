'''



'''
import bisect

def phase_mutations(transcript1, transcript2, hapdict, chromosome, start, end):
    '''
        Applies phased mutations to 2 transcripts to create the "base" for any further combinatorial mutations of sequence
        transcript1: transcript object for first chromosome
        transcript2: transcript object for second chromosome
        hapdict: haplotype dictionary
        chromosome: chromosome of transcript
        start: start of applying mutations
        end: end of applying mutations
    '''
    startindex = hapdict.bisect((chromosome, start)) - 1
    endindex = hapdict.bisect((chromosome, end)) + 1
    for x in range(startindex, endindex):
        mute_type = 'V'
        if start <= hapdict.items()[x][0][1] <= end:
            if len(hapdict.items()[x][1][3]) > len(hapdict.items()[x][1][2]):
                mute_type = 'I'
            elif len(hapdict.items()[x][1][3]) < len(hapdict.items()[x][1][2]):
                mute_type = 'D'

            if mute_type = 'V':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][3], hapdict.items()[x][0][1], mute_type)
                else:
                    transcript2.edit(hapdict.items()[x][1][3], hapdict.items()[x][0][1], mute_type)
            if mute_type = 'I':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][3][1:], hapdict.items()[x][0][1]+1, mute_type)
                else:
                    transcript2.edit(hapdict.items()[x][1][3][1:], hapdict.items()[x][0][1]+1, mute_type)
            elif mute_type = 'D':
                if hapdict.items()[x][1][0] == 1:
                    transcript1.edit(hapdict.items()[x][1][2][1:], hapdict.items()[x][0][1]+1, mute_type)
                else:
                    transcript2.edit(hapdict.items()[x][1][2][1:], hapdict.items()[x][0][1]+1, mute_type)
            else:
                print "Mutation type error in phase_mutations"
    
    return transcript1, transcript2


def combinatorics_mutations(transcript_list, exon_list, chromosome, start, end, vcf, freqpos, combinations_limit = None)
    '''
        Takes non-phased mutations from VCF and applies to all transcripts that apply 
    '''
    masterseqlist = []
    masterallelelist = []
    mastertranscripttag = []
    transcript_tag = 0
    for transcripts in transcript_list:
        sequence_list = []
        allele_list = []
        sequence_list.append(transcript.seq(start, end, genome=True))
        allele_list.append(transcript.get_freq(start, end, genome=True))
        
        vcffile = open(vcf, "r")
        for line in vcffile:
            if not line or line[0] == '#':
                continue
            stripped = line.strip().split()
            in_exon = 0
            for (exonstart, exonend) in exonlist:
                if int(stripped[1]) >= exonstart and int(stripped[1]) <= exonend:
                    in_exon_flag = 1
            if in_exon_flag != 1:
                continue
            info = stripped[8].split(':')
            allelefreqline = stripped[9].split(':')
            freqscore = float(allelefreqline[freqpos].rstrip('%'))
            mutpos = int(stripped[1]) - start
            for sequences in sequence_list:
                if stripped[0] == chromosome and int(stripped[1]) <= end and int(stripped[1]) >= start and in_exon_flag == 1:
                    numseq = len(sequence_list)
                    for x in range(0, numseq):
                        if allele_list[x][mutpos] == 0:
                            newseq = list(sequence_list[x])
                            newallele = list(allele_list[x])
                            newallele[mutpos] = freqscore
                            newseq[mutpos] = stripped[4]
                            if len(stripped[3]) > 1:
                                for x in range(0, len(stripped[3])-1):
                                    newseq[mutpos+x+1] = ""
                        sequence_list.append(newseq)
                        allele_list.append(newallele)

        vcffile.close()
        masterseqlist.append(sequence_list)
        masterallelelist.append(allelelist)
        mastertranscripttag.append(transcript_tag)
        transcript_tag += 1
    
    return













