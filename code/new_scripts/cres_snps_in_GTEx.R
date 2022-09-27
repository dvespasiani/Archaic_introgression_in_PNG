## use this script to control how many SNPs are actually identified in GTEx as eQTLs
library(data.table)
library(magrittr)
library(dplyr)
library(purrr)
library(GenomicRanges)
library(biomaRt)
library(openxlsx)
library(ggplot2)
library(ggpubr)

## settings
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
source('scripts/reusable_functions.R')
options(width = 150)


## I/O directories
annotSnpsDir <- './Chromatin_states/SNPs_chromHMM_annotated'
originalSnpsDir <- './Original_files'
outPlotDir <- './Results/Plots/'

## read gtex v8 tissue specific cis-eqtls 
ciseQTLs = list.files('../Annotation_and_other_files/GTEx/GTEx_Analysis_v8_eQTL/',full.names = T,recursive = F)%>% 
  lapply(function(x) {
    x=fread(x,sep = '\t',header = T,select = c(
    'chr','variant_pos','ref','alt', 'gene_id','gene_name','maf','qval'
    ))%>%setnames(old=c('ref','alt'),new=c('REF','ALT'))
    x<-x[qval<0.05]
  }
)
names(ciseQTLs) = gsub("\\..*","",list.files('../Annotation_and_other_files/GTEx/GTEx_Analysis_v8_eQTL/',full.names = F,recursive = F))
ciseQTLs <- Map(mutate,ciseQTLs,tissue=names(ciseQTLs))

## liftover from hg38 to hg19
library(rtracklayer)
liftoverChain <- import.chain('../Annotation_and_other_files/lifover_chains/hg38tohg19/hg38ToHg19.over.chain')

ciseQTLs_hg19 <- lapply(
  copy(ciseQTLs),function(x)x=x%>%
  makeGRangesFromDataFrame(seqnames.field = 'chr',start.field = 'variant_pos',end.field = 'variant_pos', keep.extra.columns = T)%>%
  liftOver(liftoverChain)%>%unlist() %>% as.data.table()
)%>%rbindlist()%>%dplyr::select(c(all_of(range_keys[-3]),'REF','ALT','tissue','maf'))%>%unique()

## read original SNPs
originalSnps <- list.files(originalSnpsDir,recursive = F,full.names = T,pattern='*UVBaining*')%>%
  lapply(function(y)fread(y,sep='\t',header=T,select=c(1:5))%>%setnames(old=c(1:3),new=c(range_keys))
)
names(originalSnps) = ancestry[-2]

## read filtered a/naSNPs 
annotatedSnps <- list.files(annotSnpsDir,recursive = F,full.names = T,patter='new_*') %>%
  lapply(function(y)fread(y,sep='\t',header=T,select=c(range_keys,'chrom_state','REF','ALT','MAF'))%>% unique()
)
names(annotatedSnps) = ancestry

## look overlap with ciseQTLs
## 1) for all original SNPs
originalSnps_overlapCiseQTLs <- copy(originalSnps)%>%lapply(function(x){
  x<-x[ciseQTLs_hg19,on=c(names(x)[-3]),nomatch=0]
})%>%rbindlist()%>%unique()

numbOriginalSnps <- copy(originalSnps)%>%rbindlist()%>%unique()%>%nrow()
numbOriginalSnpseQTLs <- copy(originalSnps_overlapCiseQTLs)%>%nrow()
round((numbOriginalSnpseQTLs/numbOriginalSnps)*100,2)
# [1] 2.83

## 2) for the filtered set of a/naSNPs
annotatedSnps_overlapCiseQTLs <- copy(annotatedSnps)%>%lapply(function(x){
  x<-x[ciseQTLs_hg19,on=c(range_keys[-3],'REF','ALT'),nomatch=0]%>%unique()
})

numbAnnotSnps <- copy(annotatedSnps)%>%lapply(function(x)x[,c(1:3)]%>%unique()%>%nrow())
numbAnnotSnpseQTLs <- copy(annotatedSnps_overlapCiseQTLs)%>%lapply(function(x)x[,c(1:3)]%>%unique()%>%nrow())

numbAnnotSnpseQTLs
# $denisova
# [1] 132
# $modern_humans
# [1] 688
# $neanderthal
# [1] 242

purrr::map2(numbAnnotSnpseQTLs,numbAnnotSnps,function(x,y){
  round((x/y)*100,2)
})
# $denisova
# [1] 0.09
# $modern_humans
# [1] 0.3
# $neanderthal
# [1] 0.27



# ## look at 100bp region around ciseQTLs
# sequence <- as.list(seq(0,1000,100))

# dataDistribution <- lapply(copy(sequence),function(x){
#   region = x
#   regionCiseQTLs_hg19 <- copy(ciseQTLs_hg19)[,end:=start+region][,start:=start-region]
#   setkeyv(regionCiseQTLs_hg19,c(range_keys))
  
#   annotatedSnps_overlapCiseQTLsRegion <- copy(annotatedSnps)%>%lapply(function(x){
#     setkeyv(x,range_keys)
#     overlap<-foverlaps(x,regionCiseQTLs_hg19[,c(..range_keys)],type='within')%>%na.omit()
#     overlap <- overlap[,c('start','end'):=NULL]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
#     return(overlap)
#     })
  
#   numbAnnotSnps <- copy(annotatedSnps)%>%lapply(function(x)x[,c(1:3)]%>%unique()%>%nrow())
#   numbAnnotSnpseQTLsRegion <- copy(annotatedSnps_overlapCiseQTLsRegion)%>%lapply(function(x)x[,c(1:3)]%>%unique()%>%nrow())
#   prop <- purrr::map2(numbAnnotSnpseQTLsRegion,numbAnnotSnps,function(x,y){
#     round((x/y)*100,2)
#     })
#   deniData <- data.table(ancestry = 'denisova',prop = prop$denisova,numb = numbAnnotSnpseQTLsRegion$denisova)
#   mhData <- data.table(ancestry = 'modern_humans',prop = prop$modern_humans,numb = numbAnnotSnpseQTLsRegion$modern_humans)
#   neanData <- data.table(ancestry = 'neanderthal',prop = prop$neanderthal,numb = numbAnnotSnpseQTLsRegion$neanderthal)
#   distribution = rbind(deniData,mhData,neanData)
#   return(distribution)
# })

# names(dataDistribution) = sequence 

# dataDistribution <- Map(mutate,dataDistribution,distance=names(dataDistribution))%>%rbindlist()
# dataDistribution$distance = factor(dataDistribution$distance,levels=c(as.character(sequence)))

# pdf(paste(outPlotDir,'gtex/distributionCiseQTLSnpsByDistance.pdf',sep=''),width=5,heigh=5)
# ggplot(dataDistribution,aes(x=distance, y=prop,group=ancestry,col=ancestry))+
# geom_line()+
# scale_color_manual(values=my_palette_pop)+
# ylab('Proportion of SNPs')+ xlab('Distance from cis-eQTL SNPs in bp')+
# theme_classic()+
# theme(
#   axis.text.x = element_text(angle = 60, hjust=1),
#   legend.position='bottom'
# )
# dev.off()



# # ## look at how many snp are also found in GTEx
# # allsnps <- copy(snps_chrom_states)%>%lapply(function(x)x[,c(..range_keys,'MAF')]%>%unique())

# # snps_in_gtex<-copy(allsnps)%>%lapply(function(x)x[cis_eqtl_snps_hg19,on=c(range_keys[-3]),nomatch=0])

# # ## look at overlap snps eQTL region 5kb
# # ## and perform permutation to check depletion/enrichment
# # eqtl_region <- copy(cis_eqtl_snps_hg19)[,start:=start-100][,end:=end+100][,c(..range_keys)]%>%unique()
# # setkeyv(eqtl_region,range_keys)


# # allsnps_overlap <- copy(allsnps)%>%lapply(function(x){
# #   setkeyv(x,range_keys)
# #   overlap<-foverlaps(x,eqtl_region,type='within')%>%na.omit()
# #   overlap <- overlap[,c('start','end'):=NULL]%>%unique()%>%setnames(old=c('i.start','i.end'),new=range_keys[-1])
# # })


# # # ## create random set genomic regions
# # library(regioneR)

# # counts_overlap <-list()

# # for (i in 1:2){
# #     random_snps <- randomizeRegions(makeGRangesFromDataFrame(eqtl_region),genome='hg19',allow.overlaps=F)%>%as.data.table()
# #     random_snps <- random_snps[,c('width','strand'):=NULL][!seqnames %like% '_|Y|M|X']%>%setorderv(c(range_keys),1)
# #     setkeyv(random_snps,range_keys)
# #     overlap <- foverlaps(random_snps,eqtl_region,type='within')%>%na.omit()
# #     overlap <- overlap[,c('start','end'):=NULL]%>%unique()%>%setnames(old=c('i.start','i.end'),new=range_keys[-1])
# #     counts_overlap[[i]] <-data.table(permuted_scores=nrow(overlap))
# # }
# # permuted_counts <-rbindlist(counts_overlap)
# # # random_regions <- random_regions[,c('width','strand'):=NULL][!seqnames %like% '_|Y|M|X']%>%setorderv(c(range_keys),1)



# # snps_cres_states = copy(snps_chrom_states)%>%lapply(
# #   function(x)x=x[
# #     ,MAF:=round(MAF,2)
# #     ][
# #       chrom_state %in% cres_states & MAF>= 0.2
# #       ][
# #       ,c('chrom_state'):=NULL
# #     ]%>%unique()
# # )

# # ## get eqtls overlap 
# # eqtl_overlap = function(x){
# #   overlap = copy(x)
# #   overlap = lapply(
# #   overlap,function(x)x=x[
# #     cis_eqtl_snps_hg19,on=c(columns_to_read[c(1,2,4,5)]),nomatch=0
# #     ][
# #       ,c(1:6)
# #     ]%>%unique()
# #     )
# #   overlap = Map(mutate,overlap,pop=names(overlap))
# #   return(overlap)

# # }
# # cres_snps_eqtls = eqtl_overlap(snps_cres_states)

# # numb_highfreq_eqtls = lapply(cres_snps_eqtls,function(x)nrow(x[,c(1:3)]%>%unique()))
# # numb_highfreq_cres_snps = lapply(snps_cres_states,function(x)nrow(x[,c(1:3)]%>%unique()))

# # prop_highfreq_eqtls = purrr::map2(numb_highfreq_eqtls,numb_highfreq_cres_snps,function(x,y)x/y*100)

# # ## perform same analysis but on all variants 
# # all_snps_eqtls = eqtl_overlap(snps_chrom_states)

# # numb_allsnps_eqtls = lapply(all_snps_eqtls,function(x)nrow(x[,c(1:3)]%>%unique()))
# # numb_allsnps = lapply(snps_chrom_states,function(x)nrow(x[,c(1:3)]%>%unique()))

# # prop_allsnps_eqtls = purrr::map2(numb_allsnps_eqtls,numb_allsnps,function(x,y)x/y*100)
