## use this script to annotate the SNPs with chromatin state 
library(data.table)
library(magrittr)
library(dplyr)
library(GenomicRanges)
library(R.utils)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(rGREAT)

options(width = 150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

plot_dir = './Results/Plots/QCs/'
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'

source('./scripts/reusable_functions.R')

snps_chrom_states = list.files(input_dir,recursive = F,full.names = T,pattern='new_*') %>%
  lapply(function(y)fread(y,sep='\t',header=T)[,c(..range_keys)]%>% unique()
)
names(snps_chrom_states) = ancestry

##===========================================
### compare results with Findley et al 2021
##===========================================
findleyeQTLSNPs <- fread('./knownNeanInEUR/neanSNPseQTLs_findley2021.txt',sep='\t',header=T)[
  ,seqnames:=paste('chr',chr,sep='')
  ][
    ,start:=pos
    ][
      ,c('chr','pos'):=NULL
      ][
        ,end:=start+1
        ][
          ,c(..range_keys)
]
##======================================================
## use GREAT to assign snps to genes 
## and see if there is overlap with genes already reported 
## as being target of archaic introgression in EUR
##======================================================
library(rGREAT)

getGREATTargetGenes=function(t,bk){
  test=copy(t) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  background = copy(bk)%>% makeGRangesFromDataFrame()
  results=submitGreatJob(
    gr=test,
    bg=background,
    species = "hg19",
    rule= "twoClosest",
    adv_twoDistance = 1000,
    includeCuratedRegDoms = T,
    request_interval=10
    )
  targetGenes = plotRegionGeneAssociationGraphs(results,type=1,request_interval = 10)%>%as.data.table()
  return(targetGenes)
}

findleyTargetGenes <- getGREATTargetGenes(findleyeQTLSNPs,rbind(findleyeQTLSNPs[,c(..range_keys)],rbindlist(snps_chrom_states)))

myTargetGenes <- lapply(snps_chrom_states[c(1,3)],function(x){
    x<-getGREATTargetGenes(x,rbind(findleyeQTLSNPs,x))
})

library(VennDiagram)

listOfTargetGenes <- list(
    unique(na.omit(findleyTargetGenes$gene)),
    unique(na.omit(myTargetGenes[[1]]$gene)),
    unique(na.omit(myTargetGenes[[2]]$gene))
)
names(listOfTargetGenes) = c('findleyetal2021','denisova','neanderthal')

colors = c('#e76f51',my_palette_pop[[1]],my_palette_pop[[3]])

venn.diagram(
    x = listOfTargetGenes,
    category.names = c("findley2021",'denisovan','neanderthal'),
    filename = paste(plot_dir,'vennFindleyNeanDeniTargetGenes.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = colors,  
    fill = c(alpha(colors[[1]],0.3),alpha(colors[[2]],0.3),alpha(colors[[3]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = colors
)


## calculate enrichment hypergeometric test (easy way to understand phyper: https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tothumanGenes <- length(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

overlap <- length(intersect(listOfTargetGenes[[1]],listOfTargetGenes[[3]]))
set1 <-  length(listOfTargetGenes[[1]])
set2 <-  length(listOfTargetGenes[[3]])

phyper(q = overlap-1,m=set1, n=tothumanGenes-set1, k=set2,lower.tail = F, log.p = FALSE)
# p-value < 2.2e-16

calculateEnrichment <- function(el1,el2){
 overlap <- length(intersect(el1,el2))
 set1 <-  length(el1)
 set2 <-  length(el2)
 p <- phyper(q = overlap-1,m=set1, n=tothumanGenes-set1, k=set2,lower.tail = F, log.p = FALSE)
 return(p)
}

findleyDeni <- calculateEnrichment(listOfTargetGenes[[1]],listOfTargetGenes[[2]])
# [1] 1.031851e-225
findleyNean <-calculateEnrichment(listOfTargetGenes[[1]],listOfTargetGenes[[3]])
# [1] 0



### read target genes and do the same
snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
immune_cells = c('TCells','BCells')

immune_cres_snps = read_cres(snps_input_dir)%>%lapply(
  function(x)x[allfreq!='low'][cell_line%in% immune_cells]%>%setnames(old=11,new='aoi')%>%
  dplyr::select(-c(contains('element'),contains('state'),contains('cell')))%>%unique()
)
names(immune_cres_snps) = ancestry

great_enrichment=function(test){
  test=copy(test) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  background = copy(immune_cres_snps)%>%rbindlist() %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  results=submitGreatJob(
    gr=test,
    bg=background,
    species = "hg19",
    rule= "twoClosest",
    adv_twoDistance = 1000,
    includeCuratedRegDoms = T,
    request_interval=10
    )
  go_enrichment = getEnrichmentTables(results,ontology='GO Biological Process')[[1]]%>%as.data.table()
  target_genes = plotRegionGeneAssociationGraphs(results,type=1,request_interval = 10)%>%as.data.table()
  return(list(go_enrichment,target_genes))
}

densisova_go = great_enrichment(immune_cres_snps[[1]])
neanderthal_go = great_enrichment(immune_cres_snps[[3]])

target_genes = list(densisova_go[[2]],neanderthal_go[[2]])
names(target_genes) = ancestry[-2]

## 2) plot venn with genes shared among ancestries
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='https://grch37.ensembl.org') 
get_ensembl_ids = function(x){
  df = copy(x)
  df=df$gene%>%unique()
  df=getBM(
    attributes=c('hgnc_symbol','ensembl_gene_id','go_id'),
    filters = 'hgnc_symbol', 
    values=df, 
    mart = ensembl
    )%>%as.data.table()
  return(df)
}
ensembl_ids_targetgenes = lapply(list(densisova_go[[2]],neanderthal_go[[2]]),function(y)get_ensembl_ids(y))
names(ensembl_ids_targetgenes) = c('denisovan','neanderthal')

signif_GOs = list(densisova_go[[1]],neanderthal_go[[1]])%>%lapply(function(x)x=x[Hyper_Adjp_BH<0.01])
names(signif_GOs) = ancestry[-2]
signif_GOs = Map(mutate,signif_GOs,pop=ancestry[-2])

targets_w_signf_GO=copy(ensembl_ids_targetgenes)
targets_w_signf_GO=purrr::map2(targets_w_signf_GO,signif_GOs,function(x,y)x=x[x$go_id %in% y$ID])

targets_w_signf_GO=purrr::map2(target_genes,targets_w_signf_GO,function(x,y)full_join(x,y,by=c('gene'='hgnc_symbol')))
targets_w_signf_GO= copy(targets_w_signf_GO)%>%lapply(function(x)x=x[,signif_go:=ifelse(is.na(x$ensembl_gene_id)==T,'ns','s')][,c('ensembl_gene_id','go_id'):=NULL]%>%unique())


secondListTargetGenes <- list(
  unique(na.omit(findleyTargetGenes$gene)),
   unique(na.omit(targets_w_signf_GO[[1]]$gene)),
   unique(na.omit(targets_w_signf_GO[[2]]$gene))
)
names(secondListTargetGenes) = c('findley2021','denisovan','neanderthal')

venn.diagram(
    x = secondListTargetGenes,
    category.names = c("findley2021",'denisovan','neanderthal'),
    filename = paste(plot_dir,'vennFindleyNeanDeniImmuneTargetGenes.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = colors,  
    fill = c(alpha(colors[[1]],0.3),alpha(colors[[2]],0.3),alpha(colors[[3]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = colors
)

