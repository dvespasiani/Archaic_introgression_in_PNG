
## script used to perform GO enrichment analysis using rGREAT
## on TFBS SNPs annotated in CREs with MIAF/DAF >= 0.2
library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(GenomicRanges)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(ggrepel)
library(RColorBrewer)
library(viridis)

table_dir='./Results/Tables/'
plot_dir='./Results/Plots/PNG_Windo_shared_snps/'

snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
immune_cells = c('TCells','BCells')

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
source('./scripts/reusable_functions.R')
options(width = 150)

##====================
## get immune snps
##====================
cresSnps = read_cres(snps_input_dir)%>%lapply(
  function(x)x[allfreq!='low']%>%setnames(old=11,new='aoi')%>%dplyr::select(-c(contains('element'),contains('state')))
)
names(cresSnps) = ancestry

immuneCresSnps = copy(cresSnps)%>%lapply(
  function(x)x[allfreq!='low'][cell_line%in% immune_cells][,c(..range_keys,'allfreq','MAF')]%>%unique()
)

snpsOfInterest = copy(immuneCresSnps)

## get MAF/DIAF of the alleles targeting DE genes in PNG 
## and then see how many are within W indonesia and compare the allele frequencies
denisovaInWIndo <- fread('./Original_files_per_region/Denisova/deni_w_indo.gz',sep='\t',header = T,drop=c(11:14))
neanderthalInWIndo <- fread('./Original_files_per_region/Neandertal/nean_w_indo.gz',sep='\t',header = T,drop=c(11:14))

deniaSNPsInWIndo <- copy(denisovaInWIndo)[POP_ARCH_REF+POP_ARCH_ALT>1][,c(1:5)]
neanaSNPsInWIndo <- copy(neanderthalInWIndo)[POP_ARCH_REF+POP_ARCH_ALT>1][,c(1:5)]

humanSNPsWIndo <- rbind(denisovaInWIndo,neanderthalInWIndo)[
!POP_ARCH_REF+POP_ARCH_ALT>0
  ][
    !rbind(neanaSNPsInWIndo,deniaSNPsInWIndo),on=c('CHR','FROM','TO','REF','ALT')
    ][
      ,c(1:5)
]

snpsWIndo <- list(deniaSNPsInWIndo,humanSNPsWIndo,neanaSNPsInWIndo)%>%lapply(function(x)x%>%setnames(old=c(1:3),new=c(range_keys)))
names(snpsWIndo) = ancestry


sharedSnps <- purrr::map2(snpsOfInterest,snpsWIndo,function(x,y){
  copy(x)[copy(y),on=c(range_keys),nomatch=0][
    ,type:='shared'
    ]
})

numbSharedSnps = lapply(sharedSnps,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
numbSnpsOfInterest = lapply(snpsOfInterest,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
propShared = purrr::map2(numbSharedSnps,numbSnpsOfInterest,function(x,y)x/y)%>%
lapply(function(z)data.table(fraction=z,type='shared'))

propSpecific = copy(propShared)%>%lapply(function(x)x=x[,fraction:=1-fraction][,type:='specific'])

propSpecific = Map(mutate,propSpecific,pop=ancestry)%>%rbindlist()
propShared = Map(mutate,propShared,pop=ancestry)%>%rbindlist()

propSnpsWIndoPNG <- rbind(propShared,propSpecific)

pdf(paste(plot_dir,'propSnpsOfInterest_sharedWIndoPNG.pdf',sep=''),width=7,height = 5)
ggplot(propSnpsWIndoPNG, aes(x=pop,y=fraction,fill=type)) +
geom_bar(stat='identity')+
scale_fill_manual(
  values=c('#FEDC97','#B5B682'),
  name=' ',labels=c('Shared','Specific')
  )+
  ylab('\n Fraction of variants \n')+ xlab('\n \n')+ 
  theme(
    panel.spacing=unit(1, "lines"),
    panel.background =element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(),
    legend.title = element_text(),
    legend.margin = margin(c(0.5, 2, 8, 25)),
    legend.spacing.x = unit(0.5, 'cm'),
    axis.text.y = element_text(),
    axis.title.y = element_text(hjust=0.5),
    axis.text=element_text(),
    axis.line = element_line(color = "black", size = 0.5, linetype = "solid")
  )
dev.off()

## look at distribution Fst values
fst <- list.files('./fst',recursive=F,full.names=T)%>%lapply(
  function(x)fread(x,sep='\t')%>%setnames(old=c(1,2),new=c(range_keys[-3]))
)%>%rbindlist()
fst <- fst[,seqnames:=paste('chr',seqnames,sep='')]

snpsOfInterestFst <- copy(snpsOfInterest)%>%lapply(function(x)x[fst,on=range_keys[-3],nomatch=0])
snpsOfInterestFst <- Map(mutate,snpsOfInterestFst,pop=ancestry)%>%rbindlist()

plotFst<-function(x){
  comparisons = list(
    c('denisova','modern_humans'),
    c('neanderthal','modern_humans'),
    c('denisova','neanderthal')
    )

  p <- ggplot(x, aes(x=WEIR_AND_COCKERHAM_FST,fill=pop)) +
  geom_density(alpha=0.5)+
#   geom_violin(trim=T,scale = "width")+
  # geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
  scale_fill_manual(values=my_palette_pop)+
  ylab('Density')+ xlab('Fst')+ 
  # stat_compare_means(
  #   method = "wilcox.test",
  #   comparisons = comparisons,
  #   size=5
  #   )+
  theme_classic()+
    theme(
      legend.text = element_text(),
      legend.title = element_text(),
      legend.position = 'bottom',
      legend.margin = margin(c(0.5, 2, 8, 25)),
      legend.spacing.x = unit(0.5, 'cm'),
      axis.text.y = element_text(),
      axis.title.y = element_text(hjust=0.5),
      axis.text=element_text(),
      axis.line = element_line(color = "black", size = 0.5, linetype = "solid")
    )
  return(p)
}

pdf(paste(plot_dir,'distributionFstValuesSnpsOfInterest_WIndoPNG.pdf',sep=''),width=8,height = 5)
plotFst(snpsOfInterestFst)
dev.off()

# pdf(paste(plot_dir,'distribution_FstValuesImmuneCres.pdf',sep=''),width=8,height = 5)
# plotFst(cresSnpsFst[cell_line%in% immune_cells])
# dev.off()
















# ## snps in immune-related cres and no low freq
# nolowfreq_cres_snps = read_cres(snps_input_dir)%>%lapply(
#   function(x)x[allfreq!='low']%>%setnames(old=11,new='aoi')%>%dplyr::select(-c(contains('element'),contains('state')))
# )
# names(nolowfreq_cres_snps) = ancestry

# immune_nolowfreq_cres_snps = copy(nolowfreq_cres_snps)%>%lapply(
#   function(x)x[cell_line%in% immune_cells]
# )
# immune_nolowfreq_cres_snps <- immune_nolowfreq_cres_snps[c(1,3)]
# ##------------------------------------------------------
# ## compare allele frequencies across ISEA for tfbs-snps
# ##------------------------------------------------------

# ## get MAF/DIAF of the alleles targeting DE genes in PNG 
# ## and then see how many are within W indonesia and compare the allele frequencies
# denisova_in_w_indo = fread('./Original_files_per_region/Denisova/deni_w_indo.gz',sep='\t',header = T,drop=c(11:14))
# neanderthal_in_w_indo = fread('./Original_files_per_region/Neandertal/nean_w_indo.gz',sep='\t',header = T,drop=c(11:14))

# # modern_humans_w_indo = rbind(
# #   denisova_in_w_indo[POP_ARCH_REF+POP_ARCH_ALT==0],
# #   neanderthal_in_w_indo[POP_ARCH_REF+POP_ARCH_ALT==0])%>%
# #   setorderv(c('CHR','FROM'),1)%>%unique()

# ## retrieve the continuous MIAF/DAF values (not the rounded ones)
# # original_snps = list.files('./filtered_files',full.names = T,recursive = F,pattern='new_*') %>% lapply(function(y)fread(y,sep='\t',header = T))
# # names(original_snps)=c('denisova','modern_human','neanderthal')

# original_snps = list.files('./Original_files',full.names = T,recursive = F,pattern='SNPstats*') %>% lapply(
#   function(y)fread(y,sep='\t',header = T)%>%setnames(old=c(1:3),new=c(range_keys))
# )
# names(original_snps)=c('denisova','neanderthal')

# immune_nolowfreq_cres_snps_original = purrr::map2(
#   immune_nolowfreq_cres_snps,original_snps,
#   function(x,y)
#   copy(x)[copy(y),on=c(range_keys,'REF','ALT'),nomatch=0]
# )

# ## recompute (not binned) miaf
# immune_nolowfreq_cres_snps_original <- lapply(immune_nolowfreq_cres_snps_original,
#   function(x) x<-x[
#     ,MAF := ifelse(
#       POP_ARCH_REF > POP_ARCH_ALT,
#       POP_ARCH_REF / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
#       POP_ARCH_ALT / (POP_ARCH_REF + POP_ARCH_ALT + POP_NOTARCH_REF + POP_NOTARCH_ALT)
#       )
#     ][
#       ,c(..range_keys,'MAF','aoi','REF','ALT')
#     ]%>%unique()
# )


# ## now first calculate proportion shared 
# ## NB: u need to keep Windo aSNPs with instances in archaic haps
# denisova_in_w_indo_asnps = copy(denisova_in_w_indo)[POP_ARCH_REF+POP_ARCH_ALT>0]
# neanderthal_in_w_indo_asnps = copy(neanderthal_in_w_indo)[POP_ARCH_REF+POP_ARCH_ALT>0]

# # windo_asnps_nasnps = list(denisova_in_w_indo_asnps,neanderthal_in_w_indo_asnps)%>%
# # lapply(
# #   function(x)x%>%setnames(old=c(1:3),new=c(range_keys))
# # )
# # names(windo_asnps_nasnps)=c('Denisova','Neanderthal')

# windo_asnps <- list(denisova_in_w_indo_asnps,neanderthal_in_w_indo_asnps)%>%lapply(function(x){
#   x%>%setnames(old=c(1:3),new=c(range_keys))
# })
# names(windo_asnps)=c('Denisova','Neanderthal')

# shared_variants = purrr::map2(immune_nolowfreq_cres_snps_original,windo_asnps,function(x,y)
# copy(x)[copy(y),on=c(range_keys,'REF','ALT'),nomatch=0][
#     ,type:='shared'
#     ]
# )
# # shared_variants_archaics=copy(shared_variants)%>%lapply(function(x)x=x[,allele_of_interest:=ifelse(allele_of_interest==ALT,'alt','ref')])
# # shared_variants_nonarchaics=copy(shared_variants[[2]])

# # shared_variants=list(shared_variants_archaics[[1]],shared_variants_archaics[[2]])
# # names(shared_variants)= c('Denisova','Neanderthal')

# numb_shared_snps = lapply(shared_variants,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
# numb_immunecres_snps = lapply(immune_nolowfreq_cres_snps,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

# prop_shared = purrr::map2(numb_shared_snps,numb_immunecres_snps,function(x,y)x/y)%>%
# lapply(function(z)data.table(fraction=z,type='shared'))

# prop_specific = copy(prop_shared)%>%lapply(function(x)x=x[,fraction:=1-fraction][,type:='specific'])

# prop_specific = Map(mutate,prop_specific,pop=names(prop_specific))%>%rbindlist()
# prop_shared = Map(mutate,prop_shared,pop=names(prop_shared))%>%rbindlist()

# prop_snps_windo_png = rbind(prop_shared,prop_specific)

# pdf(paste(plot_dir,'fraction_highfreq_tfbs_shared_snps_png_windo.pdf',sep=''),width=7,height = 7)
# ggplot(prop_snps_windo_png, aes(x=pop,y=fraction,fill=type)) +
# geom_bar(stat='identity')+
# scale_fill_manual(
#   values=c('#FEDC97','#B5B682'),
#   name=' ',labels=c('Shared','Specific')
#   )+
#   ylab('\n Fraction of variants \n')+ xlab('\n \n')+ 
#   theme(
#     panel.spacing=unit(1, "lines"),
#     panel.background =element_rect(fill = 'white', colour = 'black',size=1),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     legend.text = element_text(),
#     legend.title = element_text(),
#     legend.margin = margin(c(0.5, 2, 8, 25)),
#     legend.spacing.x = unit(0.5, 'cm'),
#     axis.text.y = element_text(),
#     axis.title.y = element_text(hjust=0.5),
#     axis.text=element_text(),
#     axis.line = element_line(color = "black", size = 0, linetype = "solid")
#   )
# dev.off()

# ## look tissues carrying enrichment of png-specific snps
# png_specific_variants = purrr::map2(tfbs_cres_original,windo_asnps_nasnps,function(x,y)x[!y,on=c(range_keys,'REF','ALT')][,type:='png_specific'])

# enrichment_specfic_variants = copy(png_specific_variants)%>%lapply(
#     function(y)
#     y=y[
#     ,numbsnps:= .N
#     ][
#       ,numbsnps_in_tissue := .N,by=.(cell_line)
#       ][
#         ,numbsnps_notin_tissue := numbsnps - numbsnps_in_tissue
#         ][
#           ,c('cell_line','numbsnps_in_tissue','numbsnps_notin_tissue')
#           ]%>%unique()
# )

# denisovan_odds_ratio = calculate_odds_ratio(enrichment_specfic_variants[[1]],enrichment_specfic_variants[[2]],'cell_line',T,"numbsnps_in_tissue")[,ancestry:='Denisova']%>%adjust_pvalues()
# neanderthal_odds_ratio = calculate_odds_ratio(enrichment_specfic_variants[[3]],enrichment_specfic_variants[[2]],'cell_line',T,"numbsnps_in_tissue")[,ancestry:='Neanderthal']%>%adjust_pvalues()

# asnps_odds_ratio = rbind(denisovan_odds_ratio,neanderthal_odds_ratio)

# ## plot the results
# pdf(paste(plot_dir,'enrichment_png_specific_tfbs_snps.pdf',sep=''),width=8,height = 5)
# ggplot(asnps_odds_ratio, aes(x=elements, y = odds_ratio,label = p_signif))+ 
#   geom_point(aes(colour = elements))+
#   geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
#   scale_colour_manual(values = tissue_colors)+
#   geom_hline(yintercept=1,linetype='dashed',size=.5)+
#   geom_text(aes(y= asnps_odds_ratio$upper_ci+0.05),size=5)+
#   xlab(' ')+ylab('OR aSNPs vs naSNPs')+
#   facet_wrap(~ancestry)+
#   theme(
#           panel.background =element_rect(fill = 'white', colour = 'black',size=1),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           strip.text.y = element_text(hjust = 0.5),
#           strip.background = element_rect(color = 'black', linetype = 'solid'),
#           strip.background.y = element_blank(),
#           strip.background.x =element_blank(),
#           legend.position = "bottom",
#           legend.key = element_rect(fill = "white", colour = "black"),
#           axis.line = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank()
#           )
# dev.off()

# ## then calculate the MIAF/DAF for variants shared with Windo
# shared_variants = lapply(
#   shared_variants,function(x)x=x[
#     ,Windo_MIAF:=ifelse(
#       allele_of_interest=='alt',
#       POP_ARCH_ALT / (POP_ARCH_REF + POP_ARCH_ALT + POP_NOTARCH_REF + POP_NOTARCH_ALT),
#       POP_ARCH_REF / (POP_ARCH_REF+ POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)
#       )
#       ][
#         ,Windo_DAF:=ifelse(
#           allele_of_interest=='alt',
#           POP_NOTARCH_ALT / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
#           POP_NOTARCH_REF / (`POP_ARCH_REF`+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)
#           )
#           ]%>%unique()
# )

# all_shared_snps = list(
#   shared_variants[[1]][,c(..range_keys,'MAF','Windo_MIAF')],
#   # shared_variants[[2]][,c(..range_keys,'MAF','Windo_DAF')]%>%setnames(old='Windo_DAF',new='Windo_MIAF'),
#   shared_variants[[2]][,c(..range_keys,'MAF','Windo_MIAF')]
# )
# names(all_shared_snps) = c('denisova','neanderthal')
# all_shared_snps = Map(mutate,all_shared_snps,pop=names(all_shared_snps))%>%rbindlist()


# ## plot MIAF/DAF for tfbs-disrupting SNPs
# pdf(paste(plot_dir,'png_windo_snp_freq_corr.pdf',sep=''),width=5,height = 5)
# ggscatter(
#   all_shared_snps, x = "MAF", y = "Windo_MIAF",color="pop",
#   palette = my_palette_pop,facet.by= "pop",ncol=1,xlab='Allele frequency PNG',ylab='Allele frequency W Indonesia',
#   conf.int = F,cor.coef = TRUE,cor.coeff.args = list(method = "spearman")
#   )+xlim(0,1)+ylim(0,1)+geom_abline()+
#   theme(
#     panel.background =element_rect(fill = 'white', colour = 'black',size=1),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     strip.text.y = element_text(hjust = 0.5),
#     strip.background = element_rect(color = 'black', linetype = 'solid'),
#     strip.background.y = element_blank(),
#     strip.background.x =element_blank(),
#     legend.position = "bottom",
#     legend.key = element_rect(fill = "white", colour = "black"),
#     axis.line = element_blank()
#   )
# dev.off()

##---------------------
## export supp table
##---------------------
# supp_table_or_pngspecific = make_supp_table(asnps_odds_ratio,'tissue')%>%setorderv(c('ancestry','tissue'),1)
# write.xlsx(supp_table_or_pngspecific,paste(table_dir,'Supp_Table_TFBS_OR_pngspecific_pvals.xlsx',sep=''),append=T,overwrite=T)


# ## get target genes and check OAS
# target_genes=list.files('./Motifbreak/GREAT_GO_terms/target_genes/',full.names=T,recursive=F)%>%lapply(function(x)fread(x,sep='\t',header=T))
# names(target_genes)=c('Denisova','Modern_humans','Neanderthal')

# deni_targets=copy(target_genes[[1]])
# pngspecific_deni_snps=copy(png_specific_variants[[1]])

# pngspecific_deni_snps_w_targets=pngspecific_deni_snps[deni_targets,on=range_keys,nomatch=0][,cell_line:=NULL]%>%unique()

# ## check oas snps
# oas_snps=copy(pngspecific_deni_snps_w_targets)[,c('width','strand'):=NULL][gene %like% 'OAS']%>%unique()

# oas_snps_in_windo=copy(rbindlist(windo_asnps_nasnps))[oas_snps,on=range_keys,nomatch=0]
