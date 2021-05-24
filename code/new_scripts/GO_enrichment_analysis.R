## script used to perform GO enrichment analysis using rGREAT
## on TFBS SNPs annotated in CREs with MIAF/DAF >= 0.2
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


table_dir='./Results/Tables/'
plot_dir='./Results/Plots/GREAT/'

snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
columns_to_read = c(1,7:12,4,5,15)
immune_cells = c('TCells','BCells')

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
options(width = 150)

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

## read cres snps
cres_snps = read_cres(snps_input_dir,columns_to_read,immune_cells)
names(cres_snps) = gsub("\\.txt$","",list.files(snps_input_dir,recursive = F,full.names = F))

## read tfbs snps
hocomoco_tfbs = read_tfbs('./Motifbreak/output_files','hocomoco')
jaspar_tfbs = read_tfbs('./Motifbreak/output_files','hocomoco')

## combine hocomoco and jaspar results
combined_tfbs = purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
  lapply(
    function(x)x=x[
      ,c(..range_keys,'REF','ALT','geneSymbol','providerName','alleleDiff','effect')
      ] %>% unique()
)

tfbs_cres = purrr::map2(cres_snps,combined_tfbs,function(x,y)inner_join(x,y,by=c(range_keys,'REF',"ALT")))
tfbs_cres = Map(mutate,tfbs_cres,'pop'=names(tfbs_cres))
tfbs_cres = lapply(
  tfbs_cres,function(x)x=x[
    MAF>=0.2
    ][
      ,c(..range_keys)
    ]%>%unique()
)

##-----------------
## GO enrichment
##-----------------
## GREAT enrichments
great_enrichment=function(test){
  test=copy(test) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  background = copy(tfbs_cres)%>%rbindlist() %>% makeGRangesFromDataFrame(keep.extra.columns = T)
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

densisova_go = great_enrichment(tfbs_cres[[1]])
modern_human_go = great_enrichment(tfbs_cres[[2]])
neanderthal_go = great_enrichment(tfbs_cres[[3]])

## plot distance from TSS 
## if SNP is associated with multiple genes take the closest one
target_genes = list(densisova_go[[2]],modern_human_go[[2]],neanderthal_go[[2]])%>%
lapply(
  function(x)x=x[
    ,.SD[which.max(abs(distTSS))], by=.(seqnames,start,end)
  ][
    ,log10_abs_dist:= log10(abs(distTSS))
    ]
)
names(target_genes) = c('Denisova','Modern_Humans','Neanderthal')

target_genes = Map(mutate,target_genes,pop=names(target_genes))%>%rbindlist()
  
pdf(paste(plot_dir,'snp_target_distance.pdf',sep=''),width=10,height = 7)
ggplot(target_genes, aes(x=log10_abs_dist,fill=pop)) +
    geom_density(alpha=0.5)+
    scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Modern Humans','Neanderthal'))+
    xlab('Log 10 absolute distance SNP-targets')+ylab('Density')+
    theme(panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          legend.margin = margin(c(0.5, 2, 8, 25)),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black", size = 0, linetype = "solid"))
dev.off()

## export table with all significant GO enriched terms 
## and plot the top 30 most enriched terms
signif_GOs = list(densisova_go[[1]],modern_human_go[[1]],neanderthal_go[[1]])%>%lapply(function(x)x=x[Hyper_Adjp_BH<0.01])
names(signif_GOs) = c('Denisova','Modern_Humans','Neanderthal')
signif_GOs = Map(mutate,signif_GOs,pop=names(signif_GOs))

write.xlsx(signif_GOs,paste(table_dir,'Supp_Table_GO_enriched_terms.xlsx',sep=''),append=T)

library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP")

deni_human = mgoSim(signif_GOs[[1]]$ID, signif_GOs[[2]]$ID,semData=hsGO, measure="Wang", combine='BMA')
nean_human = mgoSim(signif_GOs[[3]]$ID, signif_GOs[[2]]$ID,semData=hsGO, measure="Wang", combine='BMA')
deni_nean = mgoSim(signif_GOs[[1]]$ID, signif_GOs[[3]]$ID,semData=hsGO, measure="Wang", combine='BMA')

vector_of_scores <- c(deni_human, deni_nean,nean_human)
matrix_of_scores <- matrix(0,3,3)
rownames(matrix_of_scores) = c('Denisova','Modern Human','Neanderthal')
colnames(matrix_of_scores) = c('Denisova','Modern Human','Neanderthal')

matrix_of_scores[ col(matrix_of_scores) < row(matrix_of_scores) ] <- vector_of_scores
matrix_of_scores <- matrix_of_scores + t(matrix_of_scores)
diag(matrix_of_scores) <- 1


pdf(paste(plot_dir,'GO_semantic_similarity_score.pdf',sep=''),width=10,height = 7)
ggcorrplot::ggcorrplot(
  matrix_of_scores,
  type='lower',
  lab = TRUE,
  legend.title = "Wang similarity score",
  ggtheme =theme(
    legend.background=element_rect(),
    axis.ticks=element_blank(),
    axis.text.x=element_text(angle=0, hjust=1.10),
    axis.text.y=element_text(angle=0, vjust=0.8),
    axis.title=element_text(),
    axis.line = element_line(color = "white",size = 0, linetype = "solid"),
    panel.background = element_rect(fill = 'white', size = 0.5,colour = 'white'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    title=element_text()
    )
)
dev.off()


go_enrichment_plot = function(df,cols){
  ggplot(df[1:10,c(1,2,13,14)], aes(x=reorder(name,-log10(Hyper_Adjp_BH)), y=-log10(Hyper_Adjp_BH),fill=pop)) +
  geom_bar(stat = 'identity',position = 'dodge')+
  xlab(" ") +ylab("\n -Log10 (P) \n ") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(df$Hyper_Adjp_BH)), by = 2), 1)) +
  scale_fill_manual(name= " ",values = cols)+
  theme(
    legend.position='none',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, hjust=1.10),
    axis.text.y=element_text(angle=0, vjust=0.8),
    axis.title=element_text(),
    axis.line = element_line(color = "black",size = 0, linetype = "solid"),
    panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    title=element_text()
    ) +
  coord_flip()
  
}

pdf(paste(plot_dir,'denisova_go.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(signif_GOs[[1]],my_palette_pop[1])
dev.off()

pdf(paste(plot_dir,'modernhuman_go.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(signif_GOs[[2]],my_palette_pop[2])
dev.off()

pdf(paste(plot_dir,'neanderthal_go.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(signif_GOs[[3]],my_palette_pop[3])
dev.off()

##-------------------------------------------------------
## look at genes associated with significant GO terms
##-------------------------------------------------------
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 

get_signif_targets = function(x){
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

signif_GOs_targetgenes = split(target_genes,by='pop')%>%lapply(function(y)get_signif_targets(y))

signif_GOs_targetgenes = map2(
  signif_GOs_targetgenes,signif_GOs,
  function(x,y)x=x[
    go_id %in% y$ID
    ]%>%unique()
)

## now look at how many genes are commonly targeted by the different SNPs
genes = copy(signif_GOs_targetgenes)%>%lapply(function(x)x=x[,ensembl_gene_id]%>%na.omit()%>%unique())
venn.diagram(
    x = genes,
    category.names = c("Denisova",'Modern Humans','Neanderthal'),
    filename = paste(plot_dir,'target_genes_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = my_palette_pop,  
    fill = c(alpha(my_palette_pop[[1]],0.3), alpha(my_palette_pop[[2]],0.3),alpha(my_palette_pop[[3]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27,27, 135),
    cat.dist = c(0.055, 0.055,0.055),
    cat.fontfamily = "sans",
    cat.col = my_palette_pop
)

## get genes associated with top 40 GOs for cytoscape
topGOs = lapply(signif_GOs,function(x)x=x[1:40,])
genes_topGOs = copy(signif_GOs_targetgenes)
genes_topGOs = map2(
  genes_topGOs,topGOs,
  function(x,y)x=x[
    go_id %in% y$ID
    ][
      ,go_id:=NULL
    ]%>%unique()
)

##-----------------------------------------
## check if any of (ALL) the target genes 
## is DE across ISEA between png and Windo
##------------------------------------------
count_genes = function(x){
  numb_genes = copy(x)[,'ensembl_gene_id'] %>% unique() %>% nrow() 
  return(numb_genes)
}

## read mtw kor de genes file
mtw_kor_de = fread(
  './DE_genes/DE_genes_MTW_KOR.txt.gz',sep=' ',header = F,drop='V2',
  col.names = c("ensembl_gene_id","logFC","AveExpr", "t", "P.Value","adj.P.Val", "B"))[
    adj.P.Val<=0.01
]

target_genes_de = lapply(
  signif_GOs_targetgenes,
  function(x)x=x[
    mtw_kor_de,on='ensembl_gene_id',nomatch=0
    ][
      ,go_id:=NULL
      ]%>%unique()
)

## get proportion of target genes that are also DE
## out of all target genes associated with significant GO terms
numb_target_signifGO =lapply(signif_GOs_targetgenes,function(x)count_genes(x))
numb_target_de =lapply(target_genes_de,function(x)count_genes(x))

prop_target_de = purrr::map2(numb_target_de,numb_target_signifGO,function(x,y)x/y*100)

## test if there is a significant enrichment in the numb of DE genes targeted by deni/nean/png vs all human genes
all_human_ensembl_genes = getBM(attributes='ensembl_gene_id',mart = ensembl)

enrichment_pval=function(x,y){

  numb_targets_de = copy(x)%>%count_genes()
  numb_targets = copy(y)%>%count_genes()
  tot_numb_de_genes = copy(mtw_kor_de)%>%count_genes()  

  test=phyper(
    numb_targets_de-1,
    tot_numb_de_genes,
    nrow(all_human_ensembl_genes)-tot_numb_de_genes, ## all non de genes
    numb_targets,
    lower.tail = F
    )
  return(test)
}

enrichment_pval(target_genes_de[[1]],signif_GOs_targetgenes[[1]]) # deni
enrichment_pval(target_genes_de[[2]],signif_GOs_targetgenes[[2]]) # mh
enrichment_pval(target_genes_de[[3]],signif_GOs_targetgenes[[3]]) # nean

##-----------------------------------------------
## now compare allele frequencies across ISEA
##-----------------------------------------------

## get MAF/DIAF of the alleles targeting DE genes in PNG 
## and then see how many are within W indonesia and compare the allele frequencies
denisova_in_w_indo = fread('./Original_files_per_region/Denisova/deni_w_indo.gz',sep='\t',header = T,drop=c(11:14))
neanderthal_in_w_indo = fread('./Original_files_per_region/Neandertal/nean_w_indo.gz',sep='\t',header = T,drop=c(11:14))

modern_humans_w_indo = rbind(
  denisova_in_w_indo[POP_ARCH_REF+POP_ARCH_ALT==0],
  neanderthal_in_w_indo[POP_ARCH_REF+POP_ARCH_ALT==0])%>%
  setorderv(c('CHR','FROM'),1)%>%unique()

## retrieve the continuous MIAF/DAF values (not the rounded ones)
original_snps = list.files('./filtered_files',full.names = T,recursive = F) %>% lapply(function(y)fread(y,sep='\t',header = T))
names(original_snps)=c('denisova','modern_human','neanderthal')

tfbs_cres_original = purrr::map2(tfbs_cres,original_snps,function(x,y)x[y,on=c(range_keys),nomatch=0])
tfbs_cres_original = lapply(tfbs_cres_original,function(x)x=setnames(x,old=7,new='allele_of_interest'))


## now first calculate proportion shared 
## NB: u need to keep Windo aSNPs with instances in archaic haps
denisova_in_w_indo_asnps = copy(denisova_in_w_indo)[POP_ARCH_REF+POP_ARCH_ALT>0]
neanderthal_in_w_indo_asnps = copy(neanderthal_in_w_indo)[POP_ARCH_REF+POP_ARCH_ALT>0]

windo_asnps_nasnps =list(denisova_in_w_indo_asnps,modern_humans_w_indo,neanderthal_in_w_indo_asnps)%>%
lapply(
  function(x)x%>%setnames(old=c(1:3),new=c(range_keys))
)

## now get shared variants, then calculate and plot the % of shared variants
shared_variants = purrr::map2(tfbs_cres_original,windo_asnps_nasnps,function(x,y)x[y,on=c(range_keys,'REF','ALT'),nomatch=0])

numb_shared_snps = lapply(shared_variants,function(x)x[,c(1:3)]%>%unique()%>%nrow())
numb_snps_with_targets = lapply(tfbs_cres_original,function(x)x[,c(1:3)]%>%unique()%>%nrow())

prop_shared = purrr::map2(numb_shared_snps,numb_snps_with_targets,function(x,y)x/y)%>%
lapply(function(z)data.table(fraction=z,type='shared'))

prop_specific = copy(prop_shared)%>%lapply(function(x)x=x[,fraction:=1-fraction][,type:='specific'])

prop_specific = Map(mutate,prop_specific,pop=names(prop_specific))%>%rbindlist()
prop_shared = Map(mutate,prop_shared,pop=names(prop_shared))%>%rbindlist()

prop_snps_windo_png = rbind(prop_shared,prop_specific)

pdf(paste(plot_dir,'fraction_highfreq_tfbs_shared_snps_png_windo.pdf',sep=''),width=7,height = 7)
ggplot(prop_snps_windo_png, aes(x=pop,y=fraction,fill=type)) +
geom_bar(stat='identity')+
scale_fill_manual(
  values=c('lightskyblue2','plum2'),
  name=' ',labels=c('Shared','Specific')
  )+
  ylab('\n Fraction of variants \n')+ xlab('\n \n')+ 
  theme(
    panel.spacing=unit(1, "lines"),
    panel.background =element_rect(fill = 'white', colour = 'black',size=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(),
    legend.title = element_text(),
    legend.margin = margin(c(0.5, 2, 8, 25)),
    legend.spacing.x = unit(0.5, 'cm'),
    axis.text.y = element_text(),
    axis.title.y = element_text(hjust=0.5),
    axis.text=element_text(),
    axis.line = element_line(color = "black", size = 0, linetype = "solid")
  )
dev.off()


## calculate MIAF/DAF for variants in Windo
shared_variants = lapply(
  shared_variants,function(x)x=x[
    ,Windo_MIAF:=ifelse(
      allele_of_interest=='alt',
      POP_ARCH_ALT / (POP_ARCH_REF + POP_ARCH_ALT + POP_NOTARCH_REF + POP_NOTARCH_ALT),
      POP_ARCH_REF / (POP_ARCH_REF+ POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)
      )
      ][
        ,Windo_DAF:=ifelse(
          allele_of_interest=='alt',
          POP_NOTARCH_ALT / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
          POP_NOTARCH_REF / (`POP_ARCH_REF`+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)
          )
          ]%>%unique()
)

all_shared_snps = list(
  shared_variants[[1]][,c(..range_keys,'MAF','Windo_MIAF')],
  shared_variants[[2]][,c(..range_keys,'MAF','Windo_DAF')]%>%setnames(old='Windo_DAF',new='Windo_MIAF'),
  shared_variants[[3]][,c(..range_keys,'MAF','Windo_MIAF')]
)
names(all_shared_snps) = c('denisova','modern_humans','neanderthal')
all_shared_snps = Map(mutate,all_shared_snps,pop=names(all_shared_snps))%>%rbindlist()


pdf(paste(plot_dir,'png_windo_snp_freq_corr.pdf',sep=''),width=25,height = 12)
ggplot(all_shared_snps, aes(x=MAF,y=Windo_MIAF,group=pop)) +
geom_point(aes(color=pop),size=2)+
# geom_abline() +
scale_color_manual(name= " ",values = my_palette_pop)+
xlab('Frequency PNG')+ylab('Frequency Windonesia')+
facet_wrap(pop~.,ncol=3)+
theme(
  panel.spacing=unit(1, "lines"),
  panel.background =element_rect(fill = 'white', colour = 'black',size=1),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.text = element_text(),
  legend.title = element_text(),
  legend.margin = margin(c(0.5, 2, 8, 25)),
  legend.spacing.x = unit(0.5, 'cm'),
  axis.text.y = element_text(),
  axis.title.y = element_text(hjust=0.5),
  axis.text=element_text(),
  axis.line = element_line(color = "black", size = 0, linetype = "solid")
  )
dev.off()


## get papuan specific variants
## and potential target genes they might regulate
png_specific = purrr::map2(tfbs_cres_original,shared_variants,function(x,y)x[!y,on=c(range_keys,'REF','ALT')])

png_specific_with_targets = lapply(
  png_specific,function(x)x=x[
    target_genes,on=c(range_keys),nomatch=0
  ]
  
)

OAS_snps = copy(png_specific_with_targets[[1]])[gene %like% 'OAS']
OAS_snps_in_windo =copy(OAS_snps)%>%inner_join(denisova_in_w_indo, by=c('seqnames'='CHR','start'='FROM','end'='TO','REF','ALT'))

OAS_snps_tfs = OAS_snps[
  combined_tfbs[[1]],on=c(range_keys,'REF','ALT'),nomatch=0
  ][
    ,c('width','strand','log10_abs_dist','pop','providerName','alleleDiff','effect'):=NULL
]%>%unique()

## now read the list of dbSNPs and see if any of these have associated rsids
dbsnps = fread('../Annotation_and_other_files/human_genome/reduced_hg19_dbSNPs_150.txt.gz',sep='\t')%>%
setnames(old=c('ref','alt'),new=c('REF','ALT'))

OAS_snps = OAS_snps[dbsnps,on=c(range_keys[-3],'REF','ALT'),nomatch=0]

write.table(OAS_snps[,-c(9,10,12,13)],'./Results/Tables/OAS_snps.txt',sep='\t',quote=F,row.names=F,col.names=T)


