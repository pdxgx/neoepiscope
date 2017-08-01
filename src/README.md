NEOSCAN
-----
Neoscan is an epitope-prediction software for use in cancer immunotherapy research. It is useful in checkpoint blockade therapy
and identifying targets for vaccines. Neoscan works with both missense and indel mutations. Predictions are possible with or without a
tumor-specific BAM file.

License
-----
[MIT](http://choosealicense.com/licenses/mit/) License

Other Software
-----
[netMHC]  for predicting Binding Affinities
[hapCUT 2] (https://github.com/vibansal/HapCUT2) interpretation of max-cut algorithm for phasing.

Using Neoscan
-----
GTF_pickler.py
Inputs: gtf files, name of file to be pickled into
Outputs: pickled file
Command Line (In Bash): ```python /PATH/TO/GTF_pickler.py -g PATH/TO/GTF -d PATH/TO/PICKLE_FILE.p```

bowtie_index.py
Instruction: Save to same directory as neoscan.py; contains the bowtie query sequences

neoscan.py - MAIN FILE
Command line (In Bash): ```python /PATH/TO/neoscan.py -x /PATH/TO/BOWTIE/BASENAME -v PATH/TO/VCF -d PATH/TO/PICKLE_FILE.p > OUTPUT_FILE.txt```
Inputs: vcf file, pickled dictionaries from GTF_pickler.py, path to bowtie index basename
Outputs:
  Column 1 - wild type amino acid epitope/ kmer
  Column 2 - corresponding mutant type amino acid epitope/ kmer
  Column 3 - list of tuples in format (mutation position, corresponding vcf line)
