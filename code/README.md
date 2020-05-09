## Scripts
* All analyses used the R module: R-bundle-Bioconductor/3.7-intel-2018.u4-R-3.5.1 on Spartan.
* All SNPs input/output are in GRCh37/hg19 coordinates

### command line used:
* To annotate SNPs along the genome and to retrieve continental allele frequencies (1000 GP) we used the ensembl VEP tool specifying the following parameters

```
./vep --af --af_1kg --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --plugin CADD,[path_to]/CADD_GRCh37_1.4_whole_genome_SNVs.tsv.gz --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337
```

* To download Roadmap Epigenomics mnemonics BED files 

```
wget -r --no-parent -A '*coreMarks_mnemonics.bed.gz*' https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/ 
```

## Workflow
* Sort SNPs by ancestry and filter them in order to yield a high-confidence set of archaic variants
* Annotate SNPs across genome (VEP) and chromatin functional elements, comparing SNPs impact across cells
* Evaluate the function of regulatory variants via motifbreakR retrieving DNA motifs from Jaspar 2018, HOCOMOCO v.10 and ENCODE databases.
* Retrieve the set of regulated genes and compare their functions among ancestries
* Gene regulatory networks was generated using the stringApp v 1.3.0 in Cytoscape v 3.7.1

## Main packages used
* data.table (v 1.11.8)
* magrittr (v 1.5)
* dplyr (v 0.8.5)
* GenomicRanges (v 1.32.3)
* motifbreak (v 1.10.0)
* BSgenome.Hsapiens.UCSC.hg19 (v 1.4.0)
* rGREAT (v 1.14.0)
* ComplexHeatmap (v 1.20.0)
* annotatr (v 1.8.0)





