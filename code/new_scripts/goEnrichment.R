library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(rGREAT)
library(GenomicRanges)
library(biomaRt)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(ggrepel)
library(RColorBrewer)
library(viridis)

table_dir='./Results/Tables/'
plot_dir='./Results/Plots/GREAT/'
target_genes_dir='./Motifbreak/GREAT_GO_terms/target_genes/'

snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
columns_to_read = c(1,7:12,4,5,15)
immune_cells = c('TCells','BCells')

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
options(width = 150)

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

## read cres snps
immune_cres_snps = read_cres(snps_input_dir)%>%lapply(
  function(x)x[allfreq!='low'][cell_line%in% immune_cells]%>%setnames(old=11,new='aoi')%>%
  dplyr::select(-c(contains('element'),contains('state'),contains('cell')))%>%unique()
)
names(immune_cres_snps) = ancestry

# ## read tfbs snps
# hocomoco_tfbs = read_tfbs('./Motifbreak/output_files','hocomoco')
# jaspar_tfbs = read_tfbs('./Motifbreak/output_files','jaspar')

# ## combine hocomoco and jaspar results
# combined_tfbs = purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
#   lapply(
#     function(x)x=x[
#       ,c(..range_keys,'REF','ALT','geneSymbol','providerName','alleleDiff','effect')
#       ] %>% unique()
# )

# ## use this to calculate some proportions
# all_tfbs_cres = purrr::map2(cres_snps,combined_tfbs,function(x,y)inner_join(x,y,by=c(range_keys,'REF',"ALT")))
# all_tfbs_cres = Map(mutate,all_tfbs_cres,'pop'=names(all_tfbs_cres))
# all_tfbs_cres = lapply(
#   all_tfbs_cres,function(x)x=x[
#     MAF>=0.2
#     ]
# )

# numb_all_tfbs_cres_highfreq=lapply(all_tfbs_cres,function(x)x=x[,c(1:3)]%>%unique()%>%nrow())

# ## these snps are those active in at least 1 immune cells
# immune_tfbs_cres = purrr::map2(immune_cres_snps,combined_tfbs,function(x,y)inner_join(x,y,by=c(range_keys,'REF',"ALT")))
# immune_tfbs_cres = Map(mutate,immune_tfbs_cres,'pop'=names(immune_tfbs_cres))
# immune_tfbs_cres = lapply(
#   immune_tfbs_cres,function(x)x=x[
#     MAF>=0.2
#     ][
#       ,c(..range_keys)
#     ]%>%unique()
# )

# numb_immune_tfbs_cres_highfreq=lapply(immune_tfbs_cres,function(x)x=x[,c(1:3)]%>%unique()%>%nrow())

# prop_immune_tfbs_cres_highfreq=purrr::map2(numb_immune_tfbs_cres_highfreq,numb_all_tfbs_cres_highfreq,function(x,y)x/y*100)

##-----------------
## GO enrichment
##-----------------
## GREAT enrichments
great_enrichment=function(test,bk){
  test=copy(test) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  background = copy(bk) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
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

densisova_go = great_enrichment(immune_cres_snps[[1]],rbindlist(immune_cres_snps))
modern_human_go = great_enrichment(immune_cres_snps[[2]],rbindlist(immune_cres_snps))
neanderthal_go = great_enrichment(immune_cres_snps[[3]],rbindlist(immune_cres_snps))

## 1) plot distance from TSS 
## if SNP is associated with multiple genes take the closest one
target_genes = list(densisova_go[[2]],modern_human_go[[2]],neanderthal_go[[2]])%>%
lapply(
  function(x)x=x[
    ,.SD[which.max(abs(distTSS))], by=.(seqnames,start,end)
  ][
    ,log10_abs_dist:= log10(abs(distTSS))
    ]
)
names(target_genes) = ancestry

target_genes = Map(mutate,target_genes,pop=names(target_genes))
mapply(write.table,target_genes, file = paste(target_genes_dir,names(target_genes),'.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)


# ## 2) plot venn with genes shared among ancestries
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='https://grch37.ensembl.org') 

# get_ensembl_ids = function(x){
#   df = copy(x)
#   df=df$gene%>%unique()
#   df=getBM(
#     attributes=c('hgnc_symbol','ensembl_gene_id','go_id'),
#     filters = 'hgnc_symbol', 
#     values=df, 
#     mart = ensembl
#     )%>%as.data.table()
#   return(df)
# }

# ## if u get an error here run this httr::set_config(httr::config(ssl_verifypeer = FALSE))
# ensembl_ids_targetgenes =lapply(target_genes,function(y)get_ensembl_ids(y))

## now look at how many genes are commonly targeted by the different SNPs
genes = copy(target_genes)%>%lapply(function(x)x=x[,gene]%>%na.omit()%>%unique())

venn.diagram(
    x = genes,
    category.names = c("Denisova",'Modern humans','Neanderthal'),
    filename = paste(plot_dir,'target_genes_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = my_palette_pop,  
    fill = c(alpha(my_palette_pop[[1]],0.3),alpha(my_palette_pop[[2]],0.3),alpha(my_palette_pop[[3]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = my_palette_pop
)

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

findleyTargetGenes <- great_enrichment(
  findleyeQTLSNPs,
  rbind(findleyeQTLSNPs[,c(..range_keys)],rbindlist(immune_cres_snps)[,c(..range_keys)]))

secondTargetGenes = list(findleyTargetGenes[[2]],densisova_go[[2]],neanderthal_go[[2]])%>%
lapply(
  function(x)x=x[
    ,.SD[which.max(abs(distTSS))], by=.(seqnames,start,end)
  ][
    ,log10_abs_dist:= log10(abs(distTSS))
    ]
)
names(secondTargetGenes) = c('findleyetal2021','denisova','neanderthal')
secondGenes = copy(secondTargetGenes)%>%lapply(function(x)x=x[,gene]%>%na.omit()%>%unique())

colors = c('#e76f51',my_palette_pop[[1]],my_palette_pop[[3]])

venn.diagram(
    x = secondGenes,
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

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tothumanGenes <- length(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

# overlap <- length(intersect(secondTargetGenes[[1]],secondTargetGenes[[3]]))
# set1 <-  length(secondTargetGenes[[1]])
# set2 <-  length(secondTargetGenes[[3]])

# phyper(q = overlap-1,m=set1, n=tothumanGenes-set1, k=set2,lower.tail = F, log.p = FALSE)
# # p-value < 2.2e-16

calculateEnrichment <- function(el1,el2){
 overlap <- length(intersect(el1,el2))
 set1 <-  length(el1)
 set2 <-  length(el2)
 p <- phyper(q = overlap-1,m=set1, n=tothumanGenes-set1, k=set2,lower.tail = F, log.p = FALSE)
 return(p)
}

findleyDeni <- calculateEnrichment(secondTargetGenes[[1]],secondTargetGenes[[2]])
# [1] 1.031851e-225
findleyNean <-calculateEnrichment(secondTargetGenes[[1]],secondTargetGenes[[3]])
# [1] 0


## export table with all significant GO enriched terms 
# ## and plot the top 30 most enriched terms (10 x ancestry)
# signif_GOs = list(densisova_go[[1]],modern_human_go[[1]],neanderthal_go[[1]])%>%lapply(function(x)x=x[Hyper_Adjp_BH<0.01])
# names(signif_GOs) = ancestry
# signif_GOs = Map(mutate,signif_GOs,pop=ancestry)

# write.xlsx(signif_GOs,paste(table_dir,'Supp_Table_GO_enriched_terms.xlsx',sep=''),append=T,overwrite=T)

# ## export snps with target genes + info whether they are associated with significant GO
# targets_w_signf_GO=copy(ensembl_ids_targetgenes)
# targets_w_signf_GO=purrr::map2(targets_w_signf_GO,signif_GOs,function(x,y)x=x[x$go_id %in% y$ID])

# targets_w_signf_GO=purrr::map2(target_genes,targets_w_signf_GO,function(x,y)full_join(x,y,by=c('gene'='hgnc_symbol')))
# targets_w_signf_GO= copy(targets_w_signf_GO)%>%lapply(function(x)x=x[,signif_go:=ifelse(is.na(x$ensembl_gene_id)==T,'ns','s')][,c('ensembl_gene_id','go_id'):=NULL]%>%unique())

# mapply(write.table,targets_w_signf_GO, file = paste(target_genes_dir,names(target_genes),'.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## 3) look semantic similarity between significant GO terms
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP")

deni_human = mgoSim(signif_GOs[[1]]$ID, signif_GOs[[2]]$ID,semData=hsGO, measure="Wang", combine='BMA')
nean_human = mgoSim(signif_GOs[[3]]$ID, signif_GOs[[2]]$ID,semData=hsGO, measure="Wang", combine='BMA')
deni_nean = mgoSim(signif_GOs[[1]]$ID, signif_GOs[[2]]$ID,semData=hsGO, measure="Wang", combine='BMA')

vector_of_scores <- c(deni_human, deni_nean,nean_human)
matrix_of_scores <- matrix(0,3,3)
rownames(matrix_of_scores) = c('Denisova','Modern Human','Neanderthal')
colnames(matrix_of_scores) = c('Denisova','Modern Human','Neanderthal')

matrix_of_scores[ col(matrix_of_scores) < row(matrix_of_scores) ] <- vector_of_scores
matrix_of_scores <- matrix_of_scores + t(matrix_of_scores)
diag(matrix_of_scores) <- 1

matrix_of_scores[lower.tri(matrix_of_scores)] <- NA
melted_matrix = reshape2::melt(matrix_of_scores,na.rm=T)

# vector_of_scores <- c(deni_nean)
# matrix_of_scores <- matrix(0,2,2)
# rownames(matrix_of_scores) = c('Denisova','Neanderthal')
# colnames(matrix_of_scores) = c('Denisova','Neanderthal')

# matrix_of_scores[ col(matrix_of_scores) < row(matrix_of_scores) ] <- vector_of_scores
# matrix_of_scores <- matrix_of_scores + t(matrix_of_scores)
# diag(matrix_of_scores) <- 1

# matrix_of_scores[lower.tri(matrix_of_scores)] <- NA
# melted_matrix = reshape2::melt(matrix_of_scores,na.rm=T)


pdf(paste(plot_dir,'GO_semantic_similarity_score.pdf',sep=''),width=10,height = 7)
p <- ggplot(melted_matrix, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(
#    low=viridis(1),high=viridis(20)[20],mid=viridis(20)[10],
   low=c("#184e77","#1e6091","#1a759f"),mid=c("#168aad","#34a0a4","#52b69a","#76c893"),high=c("#99d98c","#b5e48c","#d9ed92"),
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="Wang similarity score") +
  theme_minimal()+ # minimal theme
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)
    )+
    coord_fixed()
p + geom_text(
  aes(Var2, Var1, label = value), 
  color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = 'bottom',
    legend.direction = "horizontal"
    )+
    guides(
      fill = guide_colorbar(barwidth = 7, barheight = 1,
      title.position = "top", title.hjust = 0.5)
)
dev.off()

## 4) plot a single GO plot with top 10 enriched terms
top_go = copy(signif_GOs)%>%lapply(function(x)x=x[
  1:15,c(1,2,13,14)
  ][
    ,log10_adj_p:=-log10(Hyper_Adjp_BH)]
)
top_go = Map(mutate,top_go,palette=my_palette_pop)
single_go = rbindlist(top_go)%>%setorderv('log10_adj_p',1)%>%setorderv('pop',-1)

plotGOenrichment <- function(df){
  ggplot(df, aes(x=factor(df$name,levels=df$name), y=log10_adj_p,fill=pop)) +
  geom_bar(stat = 'identity',position = 'dodge',col='black',fill=df$palette)+
  xlab(" ") +ylab("\n -Log10 (P adj.) \n ") +
  theme(
    legend.position='bottom',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, hjust=1.10),
    axis.text.y=element_text(angle=0, vjust=0.8),
    axis.title=element_text(),
    axis.line = element_line(color = "black",size = 0, linetype = "solid"),
    panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),

    title=element_text()
    ) +
  coord_flip()
  
}

# pdf(paste(plot_dir,'deniGO.pdf',sep=''),width=7,height = 7)
# plotGOenrichment(top_go[[1]])
# dev.off()
# pdf(paste(plot_dir,'modernhumanGO.pdf',sep=''),width=7,height = 7)
# plotGOenrichment(top_go[[2]])
# dev.off()
# pdf(paste(plot_dir,'neanGO.pdf',sep=''),width=7,height = 7)
# plotGOenrichment(top_go[[3]])
# dev.off()

pdf(paste(plot_dir,'singleGObarplot.pdf',sep=''),width=7,height = 7)
plotGOenrichment(single_go)
dev.off()

## oas snps
oasSNPs <- copy(targets_w_signf_GO[[1]])[gene %like% "OAS2|OAS3"]
## read deni tfbs snps
deniImmuneCresTfbs <- fread('./tmp_files/denisova_nolowfreq_allcres_tfbs.txt',sep='\t',header=T)
deniImmuneCresTfbsOas <- copy(deniImmuneCresTfbs)[oasSNPs,on=c(range_keys),nomatch=0]

## fst for these SNPs
fst <- list.files('./fst',recursive=F,full.names=T)%>%lapply(
  function(x)fread(x,sep='\t')%>%setnames(old=c(1,2),new=c(range_keys[-3]))
)%>%rbindlist()
fst <- fst[,seqnames:=paste('chr',seqnames,sep='')]

deniImmuneCresTfbsOasFst <-copy(deniImmuneCresTfbsOas)[fst,on=range_keys[-3],nomatch=0]

## geo SNPs targeting OAS gene
# OAS_snps = copy(png_specific_with_targets[[1]])[gene %like% 'OAS']
# OAS_snps_in_windo =copy(OAS_snps)%>%inner_join(denisova_in_w_indo, by=c('seqnames'='CHR','start'='FROM','end'='TO','REF','ALT'))

# OAS_snps_tfs = OAS_snps[
#   combined_tfbs[[1]],on=c(range_keys,'REF','ALT'),nomatch=0
#   ][
#     ,c('width','strand','log10_abs_dist','pop','providerName','alleleDiff','effect'):=NULL
# ]%>%unique()

# ## now read the list of dbSNPs and see if any of these have associated rsids
# dbsnps = fread('../Annotation_and_other_files/human_genome/reduced_hg19_dbSNPs_150.txt.gz',sep='\t')%>%
# setnames(old=c('ref','alt'),new=c('REF','ALT'))

# OAS_snps = OAS_snps[dbsnps,on=c(range_keys[-3],'REF','ALT'),nomatch=0]

# write.table(OAS_snps[,-c(9,10,12,13)],'./Results/Tables/OAS_snps.txt',sep='\t',quote=F,row.names=F,col.names=T)
