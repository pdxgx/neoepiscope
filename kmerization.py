#Get VCF file and analyze it
#Make the kmerization function
#Send it in to the binding affinity black box for analysis

#6/7/2017- Created file and wrote basic code
#6/8/2017- Made print function

import string, math

def evaluateCodons(snippet):
    newSnippet = snippet.replace("T", "U")
    codonTable = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"Stop", "UAG":"Stop",
    "UGU":"C", "UGC":"C", "UGA":"Stop", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    return codonTable[newSnippet]

def turnToAA(nucleotideString):
    aaString = ""
    for aa in range(len(nucleotideString)//3):
        codon = evaluateCodons(nucleotideString[3*aa:3*aa+3])
        if(codon == "Stop"):
            break
        else:
            aaString += codon
    return aaString

def myPrintFunction(kmerList):
    print("WILD TYPE" + "\t" + "MUTANT TYPE")
    for wtmtPair in kmerList:
        wt,mt = wtmtPair
        print(wt + "\t" + mt)
    return None

def kmer(nucleotideString, mutatedString = "", mutationPos = -1, mutationLength = 0):
    if(len(mutatedString) == 0):
        mutatedString = nucleotideString
    #If not given, assume that the mutation is in the middle of the nucleotide
    if(mutationPos == -1):
        mutationPos = len(nucleotideString)//2
    kmerList = list()
    AA = turnToAA(nucleotideString)
    mutatedAA = turnToAA(mutatedString)
    #Calculate which amino acid is affected
    mutatedAAPos = math.ceil(mutationPos/3)
    #Loop through window sizes
    for ksize in range(8, 12):
        #firstIndex corrects for mutations in the first half of the given string
        firstIndex = max(mutatedAAPos-ksize, 0)
        #lastIndex corrects for mutations in the latter half of the given string
        lastIndex = min(mutatedAAPos, len(mutatedAA)-ksize + 1)
        for startIndex in range(firstIndex, lastIndex):
            kmerList.append((AA[startIndex:startIndex+ksize], mutatedAA[startIndex:startIndex+ksize]))
    myPrintFunction(kmerList)
    return kmerList

kmer("ACGTGTGTCCGACCAGGTTTTAAAAAACGTGTGTCGACCAGACCAGGTTTTAAAACGTGTGTGTGT", "ACGTGTGTCCGACCAGGTTTTAAAAAACGTGTGTCGACCAGACCAGGTTTTAAAACGTGTGTGTGT")