## use this script to analyse the results of MotifbreakR
## in particular look at:
## 2. Predicted impact on PMWs
## 3. Enrichment over TF clusters (Vierstra et al. 2020)

library(data.table)
library(magrittr)
library(dplyr)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggrepel)
library(openxlsx)

options(width=150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
source('./scripts/reusable_functions.R')

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'

tfbs_input_dir = './Motifbreak/output_files'
snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'

# columns_to_read = c(1,7:12,4,5,15)
immune_cells = c('BCells','TCells')
tfbs_outdir = './tmp_files/'

##-------------
## input files
##-------------
## snps in immune-related cres and no low freq
nolowfreq_cres_snps = read_cres(snps_input_dir)%>%lapply(
  function(x)x[allfreq!='low']%>%setnames(old=11,new='aoi')%>%dplyr::select(-c(contains('element'),contains('state')))
)
names(nolowfreq_cres_snps) = ancestry

immune_nolowfreq_cres_snps = copy(nolowfreq_cres_snps)%>%lapply(
  function(x)x[cell_line%in% immune_cells]
)

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

## combine tfbs with cell-type/state info
immune_tfbs_nolowfreq_cres_snps = purrr::map2(immune_nolowfreq_cres_snps,combined_tfbs,function(x,y)x[y,on=c(range_keys,'REF',"ALT"),nomatch=0])
immune_tfbs_nolowfreq_cres_snps = Map(mutate,immune_tfbs_nolowfreq_cres_snps,'pop'=ancestry)

nolowfreq_cres_snps_tfbs = purrr::map2(nolowfreq_cres_snps,combined_tfbs,function(x,y)x[y,on=c(range_keys,'REF',"ALT"),nomatch=0])
nolowfreq_cres_snps_tfbs = Map(mutate,nolowfreq_cres_snps_tfbs,'pop'=ancestry)

## get some %
numbTFBS = lapply(nolowfreq_cres_snps_tfbs,function(x)count_snps(x))
numbCREs = lapply(nolowfreq_cres_snps,function(x)count_snps(x))

propTfbsCres = purrr::map2(numbTFBS,numbCREs,function(x,y)round(x/y*100,1))
# $denisova
# [1] 40.5

# $modern_humans
# [1] 42.3

# $neanderthal
# [1] 40.7

##-----------------------------
## calculate impact on motif
##-----------------------------
motif=fread('../Annotation_and_other_files/Motifs_clusters/motifs_clusters.txt',sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
]

## pwms lost because not reported by viestra
motifdb_pwms <- copy(nolowfreq_cres_snps_tfbs)%>%rbindlist()
motifdb_pwms <- motifdb_pwms[,database:=ifelse(providerName %like% 'H11MO','hocomoco','jaspar')][,c('providerName','database')]%>%unique()

motifdb_viestra_pwms <- copy(motifdb_pwms)%>%inner_join(motif,by=c('providerName'='motif_id'))
motifdb_viestra_pwms <- motifdb_viestra_pwms[,c('providerName','database')]

lost_pwms <- copy(motifdb_pwms)[,lost:=ifelse(providerName %in% motifdb_viestra_pwms$providerName,'no','yes')][lost=='yes'][,n:=.N,by=.(database)][,c('n','database')]%>%unique()
#     n database
# 1: 63   jaspar
# 2: 69 hocomoco

## motifbreak calculates alleleDiff as PWM score of the Ref-Alt sequences
## because we are interested either into the derived or introgressed alleles
## if the introgressed/derived allele is the alt allele keep alleleDiff as it is, otherwise reverse the sign
nolowfreq_cres_snps_tfbs_cluster <-
lapply(copy(nolowfreq_cres_snps_tfbs),
  function(x){
    x=x[
    ,c('cell_type'):=NULL
    ][
        ,alleleDiff:=ifelse(aoi=='alt',alleleDiff,-alleleDiff)
        ]%>%unique()%>%inner_join(motif,by=c('providerName'='motif_id'))
  }
)

## contingency table
contingencytable <- lapply(copy(nolowfreq_cres_snps_tfbs_cluster),function(x){
  x<- x[,c(..range_keys,'geneSymbol')]%>%unique()
  counts <- table(x$geneSymbol)
  table <- data.table(counts=counts)%>%setnames(c('geneSymbol','counts'))
  return(table)
})%>%purrr::reduce(inner_join,by='geneSymbol')%>%setnames(old=c(2:4),new=ancestry)

contingencytable <- data.frame(contingencytable)
rownames(contingencytable) = contingencytable$geneSymbol
contingencytable <- dplyr::select(contingencytable,-'geneSymbol')

## chisquared test to evaluate whether there is a significant association between
## ancestries for each cluster
chisquared_test <- function(comparison){
  chisq <- chisq.test(copy(comparison))
  pval <- chisq$p.value
  residuals <- round(chisq$residuals, 3)
  residuals <- residuals[order(residuals[,1], decreasing = TRUE),]
  chisq_results <- list(pval,residuals)
  names(chisq_results) = c('pvalue','residuals')
  return(chisq_results)
}

deni_associations <- chisquared_test(contingencytable[,c(1:2)])
nean_associations <- chisquared_test(contingencytable[,c(3,2)])

archaic_associations <- list(deni_associations,nean_associations)%>%lapply(
  function(x){
    df<-as.data.frame(copy(x$residuals))
    df<-df%>%mutate('Name'=rownames(df))%>%as.data.table()%>%setnames(old=1,new='residual')
    df <- df[,modern_humans:=NULL]
    return(df)
    }
)
archaic_associations <- Map(mutate,archaic_associations,ancestry=c('denisova','neanderthal'))%>%rbindlist()%>%setorderv('Name',1)

## plot residuals
pdf(paste(plot_dir,'archaic_tf_residual_chisquared.pdf',sep=''),width=8,height = 8)
ggplot(archaic_associations[residual>1],aes(x=residual,y=Name,fill=ancestry))+
geom_histogram(stat='identity',position=  position_dodge())+xlab('Adjusted residuals')+ ylab(' ')+
scale_fill_manual(values=my_palette_pop[-2])+
facet_wrap(ancestry~.,ncol=2)+
geom_vline(xintercept=2,linetype='dashed')+
theme_classic()+
theme(
  legend.position='bottom',
  strip.background = element_blank(),
  strip.text = element_blank(),
  axis.text.y = element_text(size=11)
  )
dev.off()


# ## permutations to identify enriched/depleted motif clusters
# set.seed(2022)
# permute <- function(tfbs){
#   list_permutations = list()
#   motifs <-unique(copy(tfbs$Name))
#   observed_value <- copy(tfbs)[,c(..range_keys,'Name')][,observed_value:=.N,by=.(Name)][,c('observed_value','Name')]%>%unique()%>%split(by='Name')

#   for(i in 1:10000){
#         permuted_score <- data.table(permutation=table(sample(motifs,nrow(tfbs),replace=T)))
#         list_permutations[[i]] = permuted_score[,permutation:=i]%>%setnames(old=c(1:2),new=c('Name','permuted_score'))
#   }
#   motif_permutation <- rbindlist(list_permutations)%>%split(by='Name')

#   enrichment = purrr::map2(motif_permutation,observed_value,function(p,o){
#     zscore = (o$observed_value-mean(p$permuted_score))/sd(p$permuted_score)
#     pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#     stat =  data.table(zscore = zscore,pval=pvalue)
#     return(stat)
#     }
#   )
#   enrichment <- Map(mutate,enrichment,tf=names(enrichment))%>%rbindlist()%>%setorderv('zscore',1)
#   enrichment <- enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]
#   return(enrichment)
# }
# deni_enriched_clusters <- permute(immune_tfbs_nolowfreq_cres_snps_cluster[[1]])[zscore>0][pval<=0.05][,ancestry:='denisova']
# nean_enriched_clusters <- permute(immune_tfbs_nolowfreq_cres_snps_cluster[[3]])[zscore>0][pval<=0.05][,ancestry:='neanderthal']

# ## intersection between enriched motif clusters
# deni_nean_common_enriched_clusters <- list(copy(deni_enriched_clusters)$tf,copy(nean_enriched_clusters)$tf)

# VennDiagram::venn.diagram(
#     x = deni_nean_common_enriched_clusters,
#     category.names = c("Denisova",'Neanderthal'),
#     filename = paste(plot_dir,'deni_nean_enriched_clusters.png',sep=''),
#     output = TRUE ,
#     imagetype="png" ,
#     height = 700 , 
#     width = 700 , 
#     resolution = 400,
#     lwd = 1,
#     col = my_palette_pop[-2],  
#     fill = c(alpha(my_palette_pop[[1]],0.3),alpha(my_palette_pop[[3]],0.3)),
#     cex = 0.5,
#     fontfamily = "sans",
#     cat.cex = 0.3,
#     cat.default.pos = "outer",
#     cat.pos = c(-27,27),
#     cat.dist = c(0.055,0.055),
#     cat.fontfamily = "sans",
#     cat.col = my_palette_pop[-2]
# )

# intersect(deni_nean_common_enriched_clusters[[1]],deni_nean_common_enriched_clusters[[2]])
# # [1] "FOX/5"       "FOX/4"       "EVI1/MECOM"  "ETS/1"       "BCL6/1"      "AHR"         "Ebox/CACCTG" "MZF1"        "ARI5B"       "CREB3/XBP1" 
# # [11] "DMRT3"  


# ## plot significantly enriched clusters
# plot_enrichments <- function(enrich_res,column,palette,xlab){
#   df <- copy(enrich_res)[,column_to_reorder:=column]
#   p <- ggplot(df,aes(x=reorder(column_to_reorder,log_zscore),y=log_zscore,fill=ancestry))+
#     geom_bar(stat='identity',position =  position_dodge())+
#     scale_fill_manual(values = palette)+
#     xlab(xlab) + ylab('log zscore')+
#     theme_classic()+
#     theme(
#       legend.position = "none"
#       # axis.text.x =element_text(angle = 75, vjust =1, hjust=1),
#       # axis.ticks.x =element_blank()
#     )+coord_flip()
#   return(p)
# }


# pdf(paste(plot_dir,'deni_enriched_tf.pdf',sep=''),width=7,height = 7)
# plot_enrichments(enrich_res=deni_enriched_clusters,column=deni_enriched_clusters$tf,palette=my_palette_pop[1],xlab='Motif Cluster')
# dev.off()

# pdf(paste(plot_dir,'nean_enriched_tf.pdf',sep=''),width=7,height = 7)
# plot_enrichments(enrich_res=nean_enriched_clusters,column=nean_enriched_clusters$tf,palette=my_palette_pop[3],xlab='Motif Cluster')
# dev.off()

plotPWM <- function(x){
  comparisons = list(
  c('denisova','modern_humans'),
  c('neanderthal','modern_humans'),
  c('neanderthal','denisova')
  )
  df <- rbindlist(copy(x))[,c(..range_keys,'Name','geneSymbol','alleleDiff','pop')]%>%unique()
  p <- ggplot(df,aes(x=pop,y=alleleDiff,fill=pop))+
      scale_fill_manual(name= " ",values = my_palette_pop)+
      xlab(' ')+ylab('Delta PWM')+
      geom_violin(trim=T,scale = "width")+
      geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
      stat_compare_means(
      method = "wilcox.test",
      method.args = list(mu=0),
      comparisons = comparisons,
      # position=position_dodge(width = 0.5),
      size=5
      )+theme_classic()+
      theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
          )
  return(p)
}
## compare delta pwms between archaic ancestries in all cells
pdf(paste(plot_dir,'tfbs_deltapwm_allCells.pdf',sep=''),width = 7,height = 7)
plotPWM(nolowfreq_cres_snps_tfbs_cluster)
dev.off()

## compare delta pwms between archaic ancestries in immune related cells
pdf(paste(plot_dir,'tfbs_deltapwm_immuneCells.pdf',sep=''),width = 7,height = 7)
plotPWM(lapply(copy(nolowfreq_cres_snps_tfbs_cluster),function(x)x[cell_line%in% immune_cells]))
dev.off()

##============================
## Findley tfbs results  
##============================

findleyResults <- fread('./knownNeanInEUR/neanTFBSresults_findley2021.txt',sep='\t',header=T,select='Motif_Name')[
  ,Motif_Name:=toupper(Motif_Name)
]%>%unique()%>%setnames('geneSymbol')%>%dplyr::pull('geneSymbol')

neanTFBS <- copy(combined_tfbs[[3]])[,geneSymbol]%>%unique()

intersect(neanTFBS,findleyResults)

venncol = c(my_palette_pop[[3]],'#e76f51')
VennDiagram::venn.diagram(
    x = list(neanTFBS,findleyResults),
    category.names = c('Neanderthal','Findley et al'),
    filename = paste(plot_dir,'nean_findley_overlap.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = venncol,  
    fill = c(alpha(venncol[1],0.3),alpha(venncol[2],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27,27),
    cat.dist = c(0.055,0.055),
    cat.fontfamily = "sans",
    cat.col = venncol
)


## numb and prop snps with DeltaPWM >/< 0
impactOnDPWM <- lapply(copy(nolowfreq_cres_snps_tfbs_cluster),function(x){
  x<-x[,c(..range_keys,'geneSymbol','alleleDiff')]%>%unique()
  x<-x[
    ,dpwm:=ifelse(alleleDiff>0,'g','l')
    ][
      ,numndpwm:=.N,by=.(dpwm)
      ][
        ,numbSnps:=.N
        ][
          ,propSnps:=round(100*(numndpwm/numbSnps),2)
          ][
            ,c('numndpwm','propSnps','dpwm')
          ]%>%unique()
})

# $denisova
#    numndpwm propSnps dpwm
# 1:     7663    56.17    l
# 2:     5979    43.83    g

# $modern_humans
#    numndpwm propSnps dpwm
# 1:    14225    58.63    l
# 2:    10037    41.37    g

# $neanderthal
#    numndpwm propSnps dpwm
# 1:     4804     56.1    l
# 2:     3759     43.9    g

## export TFBS snps with all info
tfbs_outfiles <- purrr::map2(nolowfreq_cres_snps_tfbs,nolowfreq_cres_snps_tfbs_cluster,function(x,y){
    y<-copy(y)
    x <- copy(x)[,c('pop','allfreq','cell_line','cell_type'):=NULL]%>%unique()
    x<-x[,Name:=ifelse(seqnames %in% y$seqnames & start %in% y$start,y$Name,'null')]
  }
)
tfbs_outfile_names = paste(tfbs_outdir,ancestry,'_nolowfreq_allcres_tfbs.txt',sep='')
mapply(fwrite,tfbs_outfiles, file = tfbs_outfile_names,col.names = T, row.names = F, sep = "\t", quote = F)

