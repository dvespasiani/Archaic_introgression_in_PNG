## Scripts
* All analyses used the R module: R-bundle-Bioconductor/3.7-intel-2018.u4-R-3.5.1 on Spartan.
* All SNPs input/output are in GRCh37/hg19 coordinates
* Roadmap Epigenomics mnemonics BED files used were downleaded as: 

```
wget -r --no-parent -A '*coreMarks_mnemonics.bed.gz*' https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/ 
```

## Workflow
* Sort SNPs by ancestry and filter them in order to yield a high-confidence set of archaic variants
* Annotate SNPs across genome and chromatin functional elements and assess SNPs impact across cells
* Evaluate the actual function of regulatory variants via motifbreakR retrieving DNA motifs from Jaspar 2018, HOCOMOCO v.10 and ENCODE databases.