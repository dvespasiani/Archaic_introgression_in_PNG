## use this script to investigate in which tissue aSNPs are mostly enriched
## because we are mainly interested in regulatory variation, look for CREs states
library(data.table)
library(magrittr)
library(dplyr)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(openxlsx)


setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
options(width=150)
range_keys = c('seqnames','start','end')
columns_to_read = c(range_keys,'REF','ALT','chrom_state','cell_type','cell_line','MAF')

cres_states = c('1_TssA','2_TssAFlnk','3_TxFlnk','6_EnhG','7_Enh')

## specify dirs
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
plot_dir = './Results/Plots/Chromatin_State/'


snps_chrom_states = list.files(input_dir,recursive = F,full.names = T) %>%
  lapply(function(y)fread(y,sep='\t',header=T,select = columns_to_read)%>% unique()
)
names(snps_chrom_states) = gsub("\\..*","",list.files(input_dir,recursive = F,full.names = F))

snps_cres_states = copy(snps_chrom_states)%>%lapply(
  function(x)x=x[
    ,MAF:=round(MAF,2)
    ][
      chrom_state %in% cres_states & MAF >=0.2
      ][
        ,c('chrom_state','cell_type'):=NULL
        ]%>%unique()
)


## count number of SNPs in tissue and not in tissue (to compute odds ratio later)
snps_cres_states = lapply(
  snps_cres_states,function(x)
  x=x[
    ,numbsnps:= .N
    ][
      ,numbsnps_in_tissue := .N,by=.(cell_line)
      ][
        ,numbsnps_notin_tissue := numbsnps - numbsnps_in_tissue
        ][
          ,c('cell_line','numbsnps_in_tissue','numbsnps_notin_tissue')
          ]%>%unique()
)


calculate_odds_ratio = function(asnps,nasnps,merging_keys){
    table = asnps[nasnps,on=c(merging_keys),nomatch=0]
    table[is.na(table)] = 0
    table = split(table,by=merging_keys)%>%
    lapply(function(x)x=x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%
    fisher.test()
    )
    
    fisher_test_results = copy(table)%>%
    lapply(function(x)x=data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
    ))
    fisher_test_results = Map(mutate,fisher_test_results,elements=names(fisher_test_results))%>%rbindlist()
    fisher_test_results = fisher_test_results[,significance:=ifelse(p<0.05,'*','')]

    return(fisher_test_results)
}

denisovan_odds_ratio = calculate_odds_ratio(snps_cres_states[[1]],snps_cres_states[[3]],'cell_line')[,ancestry:='Denisova']
neanderthal_odds_ratio = calculate_odds_ratio(snps_cres_states[[2]],snps_cres_states[[3]],'cell_line')[,ancestry:='Neanderthal']


asnps_odds_ratio = rbind(denisovan_odds_ratio,neanderthal_odds_ratio)

## plot the results
tissue_colors = c(
  "Adipose"='darkorange3',
  "BCells"='palegreen4',
  "Brain"='goldenrod3',
  "Digestive"="plum3",
  "Epithelial"="orange1",
  "ES_cells"='hotpink4',
  "ES_derived_cells"='dodgerblue2',
  "Heart"="palevioletred2",
  'IMR90_fetal_lung_fibroblast'='red3',
  "iPSC"='mediumpurple4',
  "Mesenchymal"='hotpink3',
  "Muscle"="indianred",
  "Myosatellite"='darkorange2',
  "Neurospheres"='lightgoldenrod2',
  "Other_cells"="grey60",
  "Smooth_muscle"="hotpink1",
  "TCells"='chartreuse4',
  "Thymus"='yellow2'
)

pdf(paste(plot_dir,'cres_tissue_enrichment.pdf',sep=''),width=8,height = 5)
ggplot(asnps_odds_ratio, aes(x=elements, y = odds_ratio,label = significance))+ 
  geom_point(aes(colour = elements))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
  scale_colour_manual(values = tissue_colors)+
  geom_hline(yintercept=1,linetype='dashed',size=.5)+
  geom_text(aes(y= asnps_odds_ratio$upper_ci+0.05),size=5)+
  xlab(' ')+ylab('odds ratio \n aSNPs vs naSNPs')+
  facet_wrap(~ancestry)+
  theme(
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          legend.position = "none",
          legend.key = element_rect(fill = "white", colour = "black"),
          axis.line = element_blank(),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
          )
dev.off()

