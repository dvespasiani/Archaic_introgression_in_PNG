## Scripts
* All analyses used the modules: r/4.0.0 and web_proxy on Spartan.
* All SNPs coordinates are in GRCh37/hg19 coordinates

### Roadmap Epigenomics mnemonics BED files downloaded via:

```
wget -r --no-parent -A '*coreMarks_mnemonics.bed.gz*' https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/ 
```
## Workflow (not super detailed)
* Sort SNPs by ancestry and filter them in order to yield a high-confidence set of archaic variants
* Annotate SNPs across genome and chromatin functional elements, comparing SNPs impact across cells
* Evaluate the function of regulatory variants via motifbreakR retrieving DNA motifs from Jaspar 2018, HOCOMOCO v11 databases. 
* Combine these SNPs with TF cluster info 
* Retrieve the set of regulated genes and compare their functions among ancestries
* Gene regulatory networks was generated using the stringApp v 1.3.0 in Cytoscape v 3.7.1






