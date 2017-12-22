Neoepiscope

Neoepiscope is a tool that takes a cancer patient's seq data (DNA and RNA from tumor and germline) and hollistically predicts a set of the most likely neoantigens for possible personalized vaccine treatments.  We found that previous methods were limiting in their breadth, often accounting for only a single, unphased SNV per epitope and not fully enumerating neoepitopes resulting from indels and start/stop codon mutations.  With this in mind, Neoepiscope emphasizes a combination of efficiency, coverage, and flexibility to provide the most accurate predictions of neoantigens resulting from a variety of mutation types and combinations.

Software prerequisites:

[HapCUT2](https://github.com/vibansal/HapCUT2)

Binding affinity predictor (either [netMHCpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) or [netMHCIIpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan))

Python packages: 

  Install via pip:
  
  intervaltree, sortedcontainers
  
  Other:
  
  bowtie_index?

Installing Neoepiscope:
INSTRUCTIONS

Instructions for running Neoepiscope:

1) If normal sample column proceeds tumor sample column in the somatic VCF, swap the columns using Neoepiscope's swap functionality for compatibility with HapCUT2:

```python neoepiscope.py swap -i [INPUT SOMATIC VCF] -o [COLUMN-SWAPPED SOMATIC VCF]```

2) If including germline variation in results, merge somatic and germline VCFs using Neoepiscope's merge functionality for combined haplotypes from HapCUT2. If columns needed to be swapped in step 1, be sure to use the swapped VCF:

```python neoepiscope.py merge -s [INPUT SOMATIC VCF] -g [INPUT GERMLINE VCF] -o [MERGED VCF]```

3) Run HapCUT2 on the merged VCF:

  a) Run HapCUT2's 'extractHAIRS' to extract haplotype relevant information from the tumor BAM file:
   
   ```extractHAIRS --indels 1 --bam [TUMOR BAM] --vcf [SOMATIC OR MERGED VCF] --out [OUTPUT FRAGMENT FILE]```
   
   b) Run HAPCUT2 to produce haplotypes:
   
    ```HAPCUT2 --fragments [FRAGMENT FILE (from step 3a)] --vcf [SOMATIC OR MERGED VCF] --output [HAPLOTYPE OUTPUT]```
    
 4) Use Neoepiscope's prep functionality to add unphased mutations to haplotype output as their own haplotypes:
 
 ```python neoepiscope.py prep -v [SOMATIC OR MERGED VCF] -c [HAPCUT2'S HAPLOTYPE OUTPUT] -o [NEW HAPLOTYPE OUTPUT]```
 
 5) Use Neoepiscope's call functionality to call neoepitopes:
 
 ```python neoepiscope.py call [OPTIONS]```
