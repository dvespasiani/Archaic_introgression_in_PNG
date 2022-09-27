library(data.table)
library(magrittr)
library(dplyr)

options(width=150)

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'
snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
target_genes_dir='./Motifbreak/GREAT_GO_terms/target_genes/'
tfbs_input_dir = './Motifbreak/output_files'
# columns_to_read = c(1,7:12,4,5,15)

##------------------------
## read and combine files
##------------------------

## this file is intense
dbsnps = fread('../Annotation_and_other_files/human_genome/reduced_hg19_dbSNPs_150.txt.gz',sep='\t')%>%
setnames(old=c('ref','alt'),new=c('REF','ALT'))

## snps-target info
# snp_target_genes=list.files(target_genes_dir,recursive=F,full.names=T)%>%lapply(function(x)fread(x,sep='\t',header=T,drop='signif_go'))

## read cres snps 
# cres_snps <- read_cres(snps_input_dir)%>%lapply(
#   function(x){
#     x <- x[,c(..range_keys,'MAF','REF','ALT')]%>%unique()
#   }
# )
# names(cres_snps) = ancestry

cres_snps = list.files('./Chromatin_states/SNPs_chromHMM_annotated/',full.names = T,recursive = F,pattern='new_*') %>% lapply(
  function(y)
  y<-fread(y,sep='\t',header = T)[chrom_state %in% cres_states][,c(..range_keys,'MAF','allfreq','REF','ALT')]%>%unique()
)
names(cres_snps)=ancestry

## tfbs disrupting snps
hocomoco_tfbs = read_tfbs(tfbs_input_dir,'hocomoco')
jaspar_tfbs = read_tfbs(tfbs_input_dir,'jaspar')

## combine hocomoco and jaspar results
combined_tfbs = purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
  lapply(
    function(x)x=x[
      ,alleleDiff:=round(alleleDiff,1)
    ][
      ,c(..range_keys,'REF','ALT','geneSymbol','providerName','alleleDiff')
      ] %>% unique()
)

combined_tfbs=lapply(
  combined_tfbs,function(x)
  x=x[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start,end,REF,ALT)]
)

## combine cres with tfbs info and then add rsid
cres_snps_w_tfbs=purrr::map2(cres_snps,combined_tfbs,function(x,y)x[y,on=c(range_keys,"REF",'ALT'),nomatch=0])
cres_snps_no_tfbs=purrr::map2(cres_snps,cres_snps_w_tfbs,function(x,y)x[!y,on=c(range_keys,"REF",'ALT')][
  ,geneSymbol:="NA"
  ][
    ,providerName:="NA"
  ][
    ,alleleDiff:='NA'
  ]
)

cres_tfbs <- purrr::map2(cres_snps_w_tfbs,cres_snps_no_tfbs,function(x,y)rbind(x,y))

cres_tfbs_with_rsid=lapply(cres_tfbs,function(x)x[dbsnps,on=c(range_keys[-3],'REF','ALT'),nomatch=0])
cres_tfbs_without_rsid=lapply(cres_tfbs,function(x)x[!dbsnps,on=c(range_keys[-3],'REF','ALT')])
cres_tfbs_without_rsid = lapply(cres_tfbs_without_rsid,function(x)x=x[,rsid:='NA'])

all_cres_tfbs=purrr::map2(cres_tfbs_with_rsid,cres_tfbs_without_rsid,function(x,y)rbind(x,y))

## add target genes info
target_genes <- list.files('./Motifbreak/GREAT_GO_terms/target_genes/',full.names=T,recursive=F)%>%lapply(
  function(x)fread(x,sep='\t',drop=c('pop','width','strand'))%>%unique()
)
names(target_genes) = ancestry


all_cres_tfbs_genes <- copy(all_cres_tfbs)
all_cres_tfbs_genes=purrr::map2(all_cres_tfbs_genes,target_genes,function(x,y)merge(x,y,by=c(range_keys),all=T))

##----------------
## export files
##----------------
mapply(write.table,all_cres_tfbs_genes, file = paste(table_dir,'h',names(all_cres_tfbs_genes),'.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

