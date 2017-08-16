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
                print "Mutation error in phase_mutations"
            
    
    return


def combinatorics_mutations()
