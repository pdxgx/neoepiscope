GTF_pickler.py
Inputs: gtf files, name of file to be pickled into
Outputs: pickled file

bowtie_index.py
Instruction: Save to same directory as neoscan.py; contains the bowtie query sequences

neoscan.py - MAIN FILE
Command line call: python /PATH/TO/neoscan.py -x /PATH/TO/BOWTIE/BASENAME -v PATH/TO/VCF -d PATH/TO/PICKLED/DICTIONARY/FILE > OUTPUT_FILE.txt
Inputs: vcf file, pickled dictionaries from GTF_pickler.py, path to bowtie index basename
Outputs:
  column 1 - wild type amino acid epitope/ kmer
  column 2 - corresponding mutant type amino acid epitope/ kmer
  column 3 - list of tuples in format (mutation position, corresponding vcf line)
