neoepiscope
-----
Neoepiscope is an neoepitope-prediction software for use in cancer immunotherapy research. Neoscan works with both missense and indel mutations, phasing somatic and germlines from the same haplotype together for more comphrensive neoepitope prediction.

License
-----
[MIT](http://choosealicense.com/licenses/mit/) License



Installing Neoepiscope
-----



Using Neoepiscope
-----

Before calling any neoepitopes, run neoepiscope in ```index``` mode to prepare dictionaries of transcript data used in neoepitope prediction:

```neoepiscope index --gtf <PATH TO GTF file> --dicts <PATH TO DIRECTORY TO HOLD PICKLED DICTIONARIES>```

To call neoepitopes from somatic mutations, ensure that data for the tumor sample in your VCF file proceeds data from a matched normal sample. If it DOES NOT, run neoepiscope in ```swap``` mode to produce a new VCF:

```neoepiscope swap --input <PATH TO INPUT VCF> --output <PATH TO WRITE SWAPPED VCF>```

If you would like to include germline variation in your neoepitope prediction, ```merge``` your somatic and germline VCFs for a sample prior to running HapCUT2:

```neoepiscope merge --germline <PATH TO GERMLINE VCF> --somatic <PATH TO SOMATIC VCF> --output <PATH TO WRITE MERGED VCF>```

Next, run HapCUT2 with your merged or somatic VCF. Before calling neoepitopes, ```prep``` your HapCUT2 output to included unphased mutations as their own haplotypes:

```neoepiscope prep --vcf <PATH TO VCF> --hapcut2-output <PATH TO HAPCUT2 OUTPUT> --output <PATH TO NEW HAPCUT OUTPUT>```

Finally, ```call`` neoepitopes:

```neoepiscope call ```
