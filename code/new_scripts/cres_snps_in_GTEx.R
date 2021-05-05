## use this script to control how many SNPs are actually identified in GTEx as eQTLs
library(data.table)
library(magrittr)
library(dplyr)
library(purrr)
library(GenomicRanges)
library(biomaRt)
library(openxlsx)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(ggpubr)


setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
options(width = 150)
range_keys = c('seqnames','start','end')
columns_to_read = c(range_keys,'REF','ALT','chrom_state','cell_type','cell_line')

cres_states = c('1_TssA','2_TssAFlnk','3_TxFlnk','6_EnhG','7_Enh')

## specify dirs
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
plot_dir = './Results/Plots/'

## read gtex v8 tissue specific cis-eqtls 
cis_eqtls = list.files('../Annotation_and_other_files/GTEx/GTEx_Analysis_v8_eQTL/',full.names = T,recursive = F)%>% 
  lapply(function(x) x=fread(x,sep = '\t',header = T,select = c(
    'chr','variant_pos','ref','alt', 'gene_id','gene_name','maf','qval'
    ))[
      qval<=0.05
      ]%>%setnames(old=c('ref','alt'),new=c('REF','ALT'))
)

names(cis_eqtls) = gsub("\\..*","",list.files('../Annotation_and_other_files/GTEx/GTEx_Analysis_v8_eQTL/',full.names = F,recursive = F))
cis_eqtls = Map(mutate,cis_eqtls,tissue=names(cis_eqtls))


## liftover from hg38 to hg19
library(rtracklayer)

liftover_chain = import.chain('../Annotation_and_other_files/lifover_chains/hg38tohg19/hg38ToHg19.over.chain')

cis_eqtls_hg19 = lapply(
  cis_eqtls,function(x)x=x%>%
  makeGRangesFromDataFrame(seqnames.field = 'chr',start.field = 'variant_pos',end.field = 'variant_pos', keep.extra.columns = T)%>%
  liftOver(liftover_chain)%>%unlist() %>% as.data.table()
)

cis_eqtl_snps_hg19 = copy(cis_eqtls_hg19)%>%rbindlist()%>%dplyr::select(c(all_of(range_keys),'REF','ALT','tissue','maf'))%>%unique()

## read snps in cres 
snps_chrom_states = list.files(input_dir,recursive = F,full.names = T) %>%
  lapply(function(y)fread(y,sep='\t',header=T,select = columns_to_read)%>% unique()
)
names(snps_chrom_states) = gsub("\\..*","",list.files(input_dir,recursive = F,full.names = F))

snps_cres_states = copy(snps_chrom_states)%>%lapply(
  function(x)x=x[
    chrom_state %in% c("4_Tx","5_TxWk")
    ][
      ,c('chrom_state','cell_type'):=NULL
    ]%>%unique()
)

cres_snps = copy(snps_cres_states)%>%lapply(function(x)x=x[,cell_line:=NULL]%>%unique())
cres_snps = Map(mutate,cres_snps,pop=names(cres_snps))

## now look at how many SNPs have been identified in GTEx
cres_snps_eqtls = lapply(
  cres_snps,function(x)x=x[
    cis_eqtl_snps_hg19,on=c(columns_to_read[c(1,2,4,5)]),nomatch=0
    ][
      ,i.end:=NULL
    ]
)

lapply(cres_snps_eqtls,function(x)nrow(x[,c(1:3)]%>%unique()))
# 72 (0.12%) denisova
# 151 (0.41) neanderthal
# 1104 (0.44) modern humans


## look at maf distribution of these variants to see if they are high-freq alleles
