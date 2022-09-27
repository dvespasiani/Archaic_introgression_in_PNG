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

columns_to_read = c(range_keys,'REF','ALT','chrom_state','cell_type','cell_line','allfreq')

## specify dirs
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
plot_dir = './Results/Plots/Chromatin_State/'
table_dir='./Results/Tables/'

snps_chrom_states = list.files(input_dir,recursive = F,full.names = T,pattern = 'new_*') %>%
  lapply(function(y)fread(y,sep='\t',header=T,select = columns_to_read)%>% unique()
)
names(snps_chrom_states) = ancestry

snps_cres_states = copy(snps_chrom_states)%>%lapply(
  function(x)x=x[
      chrom_state %in% cres_states & allfreq!='low' ## frequency filter
      ][
        ,c('chrom_state','cell_type'):=NULL
        ]%>%unique()
)

## get some %
numb_cres_snps = lapply(snps_cres_states,function(x)count_snps(x))
numb_all_snps = lapply(snps_chrom_states,function(x)count_snps(x[allfreq=='high']))

prop_highfreq_cres_highfreq_snps = purrr::map2(numb_cres_snps,numb_all_snps,function(x,y)round(x/y*100,1))
prop_highfreq_cres_highfreq_snps

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

calculate_or <- function(snps_oi,snps_noi,merging_keys){
  oi_vs_noi_or <- merge(snps_oi,snps_noi,by=c(merging_keys))%>%
  dplyr::select(c(contains('numb'),contains(all_of(merging_keys))))

  fisher_test <- copy(oi_vs_noi_or)%>%split(by=c(merging_keys))%>%
  lapply(
    function(x){
    x <- x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%fisher.test()
    x <- data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
        )
    }
  )
  or_results <- Map(mutate,fisher_test,elements=names(fisher_test))%>%rbindlist()
  return(or_results)
}

denisovan_odds_ratio = calculate_or(snps_cres_states[[1]],snps_cres_states[[2]],'cell_line')[,ancestry:='Denisova']%>%adjust_pvalues()
neanderthal_odds_ratio = calculate_or(snps_cres_states[[3]],snps_cres_states[[2]],'cell_line')[,ancestry:='Neanderthal']%>%adjust_pvalues()

asnps_odds_ratio = rbind(denisovan_odds_ratio,neanderthal_odds_ratio)

## plot the results
pdf(paste(plot_dir,'cres_tissue_enrichment.pdf',sep=''),width=8,height = 5)
ggplot(asnps_odds_ratio, aes(x=elements, y = odds_ratio,label = p_signif))+ 
  geom_point(aes(colour = elements))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
  scale_colour_manual(values = tissue_colors)+
  geom_hline(yintercept=1,linetype='dashed',size=.5)+
  geom_text(aes(y= asnps_odds_ratio$upper_ci+0.05),size=5)+
  xlab(' ')+ylab('OR aSNPs vs naSNPs')+
  facet_wrap(~ancestry)+
  theme(
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
          )
dev.off()

meanNumbaSNPS <- copy(snps_chrom_states[c(1,3)])%>%lapply(function(x){
  x<-copy(x)[
      chrom_state %in% cres_states & allfreq!='low' 
      ][
    ,c(..range_keys,'cell_line','cell_type')]%>%
    unique()

  x<-x[
      ,numbsnps_celltype:=.N,by=.(cell_type)
      ][,meanNumbSnps:=round(mean(numbsnps_celltype),1),by=.(cell_line)
  ][
    ,c('cell_line','meanNumbSnps')
    ]%>%unique()
})

meanNumbaSNPS <- Map(mutate,meanNumbaSNPS,ancestry=c('Denisova','Neanderthal'))%>%rbindlist()%>%setnames(old='cell_line',new='tissue')

pdf(paste(plot_dir,'cres_tissue_numbaSNPs.pdf',sep=''),width=8,height = 3)
ggplot(meanNumbaSNPS, aes(x=tissue, y=log10(meanNumbSnps),fill=tissue)) + 
  geom_bar(stat='identity')+
  scale_fill_manual(values = tissue_colors)+
  xlab('')+ylab('log10 mean number aSNPs')+
  facet_wrap(~ancestry)+
  theme_classic()+
  theme(
        panel.background =element_rect(fill = 'white', colour = 'black',size=1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(hjust = 0.5),
        strip.background.y = element_blank(),
        strip.background.x =element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
)
dev.off()


## write spreadsheet with stat results

numbaSNPS <- copy(snps_chrom_states[c(1,3)])%>%lapply(function(x){
  x<-copy(x)[
      chrom_state %in% cres_states & allfreq!='low' 
      ][
    ,c(..range_keys,'cell_line')
    ]%>%unique()
  x<-x[
      ,numb_snps:=.N,by=.(cell_line)
      ][
    ,c('cell_line','numb_snps')
    ]%>%unique()
})

numbaSNPS <- Map(mutate,numbaSNPS,ancestry=c('Denisova','Neanderthal'))%>%rbindlist()%>%setnames(old='cell_line',new='tissue')

supp_table <- copy(asnps_odds_ratio)%>%setnames(old='elements',new='tissue')
supp_table <- supp_table[numbaSNPS,on=c('tissue','ancestry'),nomatch=0]

write.xlsx(supp_table,paste(table_dir,'Supp_Table_OR_CREs_highfreq_pvals.xlsx',sep=''),append=T,overwrite=T)

