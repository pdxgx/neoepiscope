neoepiscope
-----
`neoepiscope` is software for predicting neoepitopes from DNA sequencing (DNA-seq) data. Where most neoepitope prediction software confines attention to neoepitopes arising from at most one somatic mutation, often just an SNV, `neoepiscope` uses assembled haplotype output of [HapCUT2](https://github.com/vibansal/HapCUT2) to also enumerate neoepitopes arising from more than one somatic mutation. `neoepiscope` also takes into account frameshifting from indels and permits personalizing the reference transcriptome using germline variants.

License
-----
`neoepiscope` is licensed under the [MIT](http://choosealicense.com/licenses/mit/) license. See `LICENSE`[LICENSE] for more details.

Installing Neoepiscope
-----

First clone this repo to `/path/to/neoepiscope/repo`; then run

```pip install /path/to/neoepiscope/repo```.

To download compatible reference annotation files (hg19 and/or GRCh38) and link installations of relevant optional softwares to `neoepiscope` (e.g. netMHCpan), you will need to use our download functionality. From within `/path/to/neoepiscope/repo` run

```python setup.py download```

and respond to the prompts as relevant for your needs.

To make sure that the software is running properly, test it by running

```python setup.py test```

Using Neoepiscope
-----

If you aren't using hg19 or GRCh38 reference builds from our download functionality, you will need to download and prepare your own annotation files. Before calling any neoepitopes, run neoepiscope in ```index``` mode to prepare dictionaries of transcript data used in neoepitope prediction:

```neoepiscope index -g <GTF> -d <DIRECTORY TO HOLD PICKLED DICTIONARIES>```

Options:

```-g, --gtf```     path to GTF file

```-d, --dicts```   path to write pickled dictionaries


To call neoepitopes from somatic mutations, ensure that data for the tumor sample in your VCF file proceeds data from a matched normal sample. If it DOES NOT, run neoepiscope in ```swap``` mode to produce a new VCF:

```neoepiscope swap -i <INPUT VCF> -o <SWAPPED VCF>```

Options:

```-i, --input```   path to input VCF

```-o, --output```  path to swapped VCF


If you would like to include germline variation in your neoepitope prediction, ```merge``` your somatic and germline VCFs for a sample prior to running HapCUT2:

```neoepiscope merge -g <GERMLINE VCF> -s <SOMATIC VCF> -o <MERGED VCF>```

Options:

```-g, --germline```  path to germline VCF

```-s, --somatic```   path to somatic VCF

```-o, --output```    path to write merged VCF


Next, run HapCUT2 with your merged or somatic VCF. Before calling neoepitopes, ```prep``` your HapCUT2 output to included unphased mutations as their own haplotypes:

```neoepiscope prep -v <VCF> -c <HAPCUT2 OUTPUT> -o <ADJUSTED HAPCUT OUTPUT>```

Options:

```-v, --vcf```               path to VCF file used to generate HapCUT2 output

```-c, --hapcut2-output```    path to original HapCUT2 output

```-o, --output```            path to output file



Finally, ```call``` neoepitopes:

```neoepiscope call -x <BOWTIE INDEX> -v <VCF> -d <DICTIONARIES> -c <HAPCUT2 OUTPUT> -o <OUTPUT> [options]```

Options:

```-x, --bowtie-index```              path to bowtie index of reference genome

```-c, --merged-hapcut2-output```     path to HapCUT2 output

```-v, --vcf```                       path to VCF file used to generate HapCUT2 output

```-d, --dicts```                     path to directory containing pickled dictionaries generated in ```index``` mode

```-o, --output_file```               path to output file

```-k, --kmer-size```                 kmer size for neoepitope prediction (default 8-11 amino acids)

```-p, --affinity-predictor```        software to use for MHC binding predictions

```-a, --alleles```                   alleles to use for MHC binding predictions

```-u, --upstream_atgs```             handling of translation from upstream start codons - ("novel" (default) only, "all", "none", "reference" only)

