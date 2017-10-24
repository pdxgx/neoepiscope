Neoscan

Neoscan is a tool that takes a cancer patient's seq data (DNA and RNA from tumor and germline) and hollistically predits a set of the most likely neoantigens for possible personalized vaccine treatments.  We found that previous methods were limiting in their breadth, only accounting for a single SNV per epitope and disregarding indels and stop codon mutations entirely.  With this in mind, Neoscan emphasizes a combination of efficiency, coverage, and flexibility to provide the most accurate predictions of neoantigens.




Austin is finally here.

For calling the function:
python PATH/TO/FUNCTION -v PATH/TO/VCF/FILE -x PATH/TO/BOWTIE/INDEX -g PATH/TO/GTF/FILE -t PATH/TO/WRITABLE/TEXT/FILE


Command to swap columns in the VCF if necessary:
awk 'substr($1, 1, 2) == "##" {print $0; next} {for(i=1; i<=NF; i+=1) {if (i == NF-1) {printf $NF "\t" } else if (i == NF) {printf $(NF-1) "\n"} else {printf $i "\t"}}}'