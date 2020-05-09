## Scripts
* All analyses used the R module: R-bundle-Bioconductor/3.7-intel-2018.u4-R-3.5.1 on Spartan.

* To annotate SNPs along the genome and to retrieve continental allele frequencies (1000 GP) we used the Ensembl VEP tool specifying the following parameters: 

```
./vep --af --af_1kg --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --plugin CADD,[path_to]/CADD_GRCh37_1.4_whole_genome_SNVs.tsv.gz --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337
```

* All SNPs input/output are in GRCh37/hg19 coordinates
* Roadmap Epigenomics mnemonics BED files used were downleaded as: 

```
wget -r --no-parent -A '*coreMarks_mnemonics.bed.gz*' https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/ 
```

## Workflow
* Sort SNPs by ancestry and filter them in order to yield a high-confidence set of archaic variants
* Annotate SNPs across genome and chromatin functional elements and assess SNPs impact across cells
* Evaluate the actual function of regulatory variants via motifbreakR retrieving DNA motifs from Jaspar 2018, HOCOMOCO v.10 and ENCODE databases.