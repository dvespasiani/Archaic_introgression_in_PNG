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

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

columns_to_read = c(range_keys,'REF','ALT','chrom_state','cell_type','cell_line','MAF')

## specify dirs
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
plot_dir = './Results/Plots/Chromatin_State/'
table_dir='./Results/Tables/'


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

## get some %
numb_cres_snps = lapply(snps_cres_states,function(x)count_snps(x))
numb_all_snps = lapply(snps_chrom_states,function(x)count_snps(x[,MAF:=round(MAF,2)][MAF>=0.2]))

prop_highfreq_cres_highfreq_snps = purrr::map2(numb_cres_snps,numb_all_snps,function(x,y)round(x/y*100,1))


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

denisovan_odds_ratio = calculate_odds_ratio(snps_cres_states[[1]],snps_cres_states[[3]],'cell_line',T)[,ancestry:='Denisova']%>%adjust_pvalues()
neanderthal_odds_ratio = calculate_odds_ratio(snps_cres_states[[2]],snps_cres_states[[3]],'cell_line',T)[,ancestry:='Neanderthal']%>%adjust_pvalues()

asnps_odds_ratio = rbind(denisovan_odds_ratio,neanderthal_odds_ratio)

## plot the results
pdf(paste(plot_dir,'cres_tissue_enrichment.pdf',sep=''),width=8,height = 5)
ggplot(asnps_odds_ratio, aes(x=elements, y = odds_ratio,label = p_signif))+ 
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


## write spreadsheet with stat results
table_OR_tissues_CREs = copy(asnps_odds_ratio)[
  ,c(5,2:4,1,8,9,7)
]%>%setnames(old='elements',new='tissue')%>%setorderv('tissue',1)

write.xlsx(table_OR_tissues_CREs,paste(table_dir,'Supp_Table_OR_CREs_highfreq_pvals.xlsx',sep=''),append=T)

