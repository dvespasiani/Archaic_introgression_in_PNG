## use this script to analyse the results of MotifbreakR
## in particular look at:
## 1. Ts:Tv ratio across each tisssue
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


options(width=150)

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'

tfbs_input_dir = './Motifbreak/output_files'
snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'

range_keys=c('seqnames','start','end')

columns_to_read = c(1,7:12,4,5,15)

immune_cells = c('BCells','TCells')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))


## read snps with tfbs predictions
hocomoco_tfbs = read_tfbs(tfbs_input_dir,'hocomoco')
jaspar_tfbs = read_tfbs(tfbs_input_dir,'jaspar')

## combine hocomoco and jaspar results
combined_tfbs = purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
  lapply(
    function(x)x=x[
      ,c(..range_keys,'REF','ALT','geneSymbol','providerName','alleleDiff','effect')
      ] %>% unique()
)
## then remove duplicates and if same tfbs is disrupted across the two databases keep the one with strongest impact
combined_tfbs = lapply(combined_tfbs,function(x)x=x[
  , .SD[which.max(abs(alleleDiff))], by=.(seqnames,start,end,REF,ALT)
  ]
)

## read cres snps
cres_snps = read_cres(snps_input_dir,columns_to_read)
names(cres_snps) = gsub("\\.txt$","",list.files(snps_input_dir,recursive = F,full.names = F))

# lapply(cres_snps,function(x)x[MAF>=0.2][,c(1:3)]%>%unique()%>%nrow())
## combine tfbs with cell-type/state info

tfbs_cres = purrr::map2(cres_snps[c(1,3)],combined_tfbs,inner_join,by=c(range_keys,'REF',"ALT"))
tfbs_cres = Map(mutate,tfbs_cres,'pop'=names(tfbs_cres))

lapply(tfbs_cres,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

## retain only those with MAF >= 0.2
highfreq_tfbs_cres = copy(tfbs_cres)%>%lapply(function(x)x[MAF>=0.2])
lapply(highfreq_tfbs_cres,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

##-------------------------
## calculate Ts:Tv ratio
##-------------------------

tvts_ratio = function(x){
  y=copy(x)[
    ,dna_ref_base:=ifelse(REF=='A' |REF=='G','purine','pyrimidine')
    ][
      ,dna_alt_base:=ifelse(ALT=='A' |ALT=='G','purine','pyrimidine')
      ][
        ,tv_ts:=ifelse(dna_ref_base=='purine' & dna_alt_base=='pyrimidine','tv',
                       ifelse(dna_ref_base=='pyrimidine' & dna_alt_base=='purine','tv','ts'))
        ]
  y=y[
    ,c('seqnames','start','end','REF','ALT','tv_ts','cell_type','cell_line')
    ][
      ,numb_ts:=sum(tv_ts =='ts'),by=.(cell_type)
      ][
        ,numb_tv:=sum(tv_ts =='tv'),by=.(cell_type)
        ][
          ,ts_tv_ratio:=numb_ts/numb_tv
          ]%>%unique()
  y=y[
    ,c('tv_ts','ts_tv_ratio','cell_type','cell_line')
    ]%>% unique() 

  avg_tstv_ratio = copy(y)[
    ,c('ts_tv_ratio','cell_line')
    ]%>%unique()%>%mutate('mean_tstv'= round(mean(ts_tv_ratio),2))%>%dplyr::pull('mean_tstv')%>%unique()
   y=y[
     ,avg_tstv:=avg_tstv_ratio
   ]
  return(y)
}


all_cres_snps_tstv = lapply(cres_snps,function(x)x=tvts_ratio(x)) ## control for ts:tv
highfreq_tfbs_snps_tstv = lapply(highfreq_tfbs_cres,function(x)x=tvts_ratio(x))


## Make Ts:Tv ratio for all cells and report significance using one sample wilcoxon test
## from this remove the cell lines with less than 3 cell types
adjust_pvalues=function(x){
  df=copy(x)
  pvals_df=copy(df)
  pvals_df=pvals_df$p
  pvals_df_adjusted=p.adjust(pvals_df,'fdr')%>%as.data.table() %>% setnames('p.adj')
  pvals_df_adjusted=pvals_df_adjusted[
    ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                             ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                    ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                           ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
    ]
   df_final=cbind(df,pvals_df_adjusted)
  return(df_final)
}

## test whether Ts:Tv ratio is significantly different from 2 (genome-wide and all cres snps expectation) 
stat_test = copy(highfreq_tfbs_snps_tstv)%>%
lapply(
  function(x)x=x[
    ,numb_celltypes:=.N,by=.(cell_line)
    ][
      numb_celltypes >= 3
      ]%>%split(by='cell_line') %>%
  lapply(
    function(y)y=y%>%
    mutate('p'=wilcox.test(y$ts_tv_ratio,mu=2,exact = F,paired = F)$p.val)%>%
    mutate('statistic'=wilcox.test(y$ts_tv_ratio,mu=2,exact = F,paired = F)$stat)
  )%>%rbindlist()
) %>% lapply(function(x)adjust_pvalues(x))

## make Supp Table with pvalues
supp_table_tstv = copy(stat_test)%>%
lapply(
  function(x)x=x[
    ,c('cell_line','p','statistic','p.adj','p.signif')
  ][
    ,cell_line:=ifelse(cell_line=='iPSC','IPSC',cell_line)
  ]%>%unique()
)
supp_table_tstv = Map(mutate,supp_table_tstv,pop=names(supp_table_tstv))%>%rbindlist()%>%setorderv('cell_line',1)


## plot Ts:Tv ratios
tstv_distribution_celltypes = copy(stat_test)%>%
lapply(
  function(x)x=x[
    ,c('ts_tv_ratio','cell_type','cell_line','p.signif')
  ][
    ,cell_line:=ifelse(cell_line=='iPSC','IPSC',cell_line)
  ]%>%unique() 
)
tstv_distribution_celltypes = Map(mutate,tstv_distribution_celltypes,pop=names(tstv_distribution_celltypes))%>%rbindlist()

# denisova_tstv_immune=copy(ts_tv_tables[[1]])[pop=='denisova'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique()
# nean_tstv_immune=copy(ts_tv_tables[[1]])[pop=='neandertal'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique() 
# png_tstv_immune=copy(ts_tv_tables[[1]])[pop=='png'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique() 


# cells_tstv_ratio$cell_line=factor(cells_tstv_ratio$cell_line,
# levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
#          'iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus'))

plot_tfbs = function(x,cells,y_scale,column){
  df=copy(x)[
    ,column_to_plot:= column
    ][
      cell_line %in% c(cells)
      ]
  # df$cell_line=factor(df$cell_line,levels=cell_line_levels)
  pvals=copy(df)[,c('pop','cell_line','p.signif')] %>% unique()
  # stat_test$cell_line=factor(stat_test$cell_line,levels=cell_line_levels)
  # modern_humans = copy(df)[pop %like% 'modern']%>%group_by(cell_line) %>%summarize(int = median(column))
  
  my_palette_pop=c(
    '#C99E10', # denisova
    # '#1E656D', # modern humans
    '#9B4F0F'  # neandertal
  )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neanderthal'))
  
  ggplot(df,aes(x=pop,y=column_to_plot,fill=pop))+
    colScale+
    xlab(' ')+ylab(' ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_jitter(position=position_jitter(0.1),size=0.2)+
    geom_text(data=pvals, aes(x=pop, y=3.05, label=p.signif), col='black', size=7)+
    geom_hline(yintercept=2, linetype="dashed", color = "darkgrey",size=0.2)+
    facet_wrap(cell_line~., ncol = 4,scales = y_scale)+
    # geom_hline(data = modern_humans, aes(yintercept=int), linetype="dashed", color = "lightgrey",size=0.2)+
    theme(strip.text.x = element_text(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          legend.margin = margin(c(0.5, 2, 8, 25)),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black", size = 0, linetype = "solid"))
}


pdf(paste(plot_dir,'tv_ts_ratio_all_cells.pdf',sep=''),width = 20,height = 20)
plot_tfbs(tstv_distribution_celltypes,tstv_distribution_celltypes$cell_line,'free_y',tstv_distribution_celltypes$ts_tv_ratio)
dev.off()

pdf(paste(plot_dir,'tv_ts_ratio_imunue_cells.pdf',sep=''),width = 8,height = 4)
plot_tfbs(tstv_distribution_celltypes,immune_cells,'fixed',tstv_distribution_celltypes$ts_tv_ratio)
dev.off()

##-----------------------------
## calculate impact on motif
##-----------------------------
## motifbreak calculates alleleDiff as PWM score of the Ref-Alt sequences
## because we are interested either into the derived or introgressed alleles
## if the introgressed/derived allele is the alt allele keep alleleDiff as it is, otherwise reverse the sign
tfbs_cres_impact = copy(tfbs_cres)%>%
lapply(
  function(x)x=x[
    ,alleleDiff:=ifelse(allele_of_interest=='alt',alleleDiff,-alleleDiff)
  ]
)%>%rbindlist()
tfbs_cres_impact = tfbs_cres_impact[MAF>=0.2]

pdf(paste(plot_dir,'TFBS_effect_by_tissue.pdf',sep=''),width = 20,height = 20)
plot_tfbs(tfbs_cres_impact,tfbs_cres_impact$cell_line,'fixed',tfbs_cres_impact$alleleDiff)
dev.off()


##----------------------------------------------------------
## look enrichment across tf clusters (Vierstra et al 2020)
##----------------------------------------------------------
## first look enrichment for all TFBS-SNPs and then only for those with MIAF/DAF >= 0.2

motif=fread('../Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
]
## count snps in clusters
## NB: not all motifs get assigned because some might differ between pwm versions
snps_in_cluster = function(x){
  df = copy(x)%>%
    inner_join(motif,by=c('providerName'='motif_id'))%>%
      dplyr::select(c(all_of(range_keys),'Name','Cluster'))
  df = df[
    ,numbtfbs_snps:= .N
    ][
      ,numbtfbs_snps_in_clust:=.N,by=.(Name)
      ][
        ,numbtfbs_snps_notin_clust:= numbtfbs_snps-numbtfbs_snps_in_clust
        ][
          ,numbtfbs_snps:=NULL
        ][
          ,c(range_keys,'Cluster'):=NULL
        ]%>%unique()
  return(df)
}


tfbs_cluster = copy(tfbs_cres)%>%lapply(function(x)snps_in_cluster(x))
highfreq_tfbs_cluster = copy(tfbs_cres)%>%lapply(function(x)x=x[MAF>=0.2]%>%snps_in_cluster())

calculate_odds_ratio = function(asnps,nasnps,merging_keys,fulljoin){
    table = merge(asnps,nasnps,by=c(merging_keys),all=fulljoin)
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

## because some clusters might not be present in one of the 2 datasets and u are doing a full join
## some odds ratios might have the upper/lower ci = Inf
## set these to 0
denisova_clust_OR = calculate_odds_ratio(tfbs_cluster[[1]],tfbs_cluster[[2]],'Name',TRUE)[,ancestry:='Denisova']
neanderthal_clust_OR = calculate_odds_ratio(tfbs_cluster[[2]],tfbs_cluster[[1]],'Name',TRUE)[,ancestry:='Neanderthal']

denisova_highfreq_clust_OR = calculate_odds_ratio(highfreq_tfbs_cluster[[1]],highfreq_tfbs_cluster[[2]],'Name',TRUE)[,ancestry:='Denisova']
neanderthal_highfreq_clust_OR = calculate_odds_ratio(highfreq_tfbs_cluster[[2]],highfreq_tfbs_cluster[[1]],'Name',TRUE)[,ancestry:='Neanderthal']


## now prepare data.table for plotting
reshape_OR = function(x,tfbs_cluster){
  df = copy(x)%>% mutate_if(is.numeric, function(y) ifelse(is.infinite(y), 1, y))%>%as.data.table()
  df = df[
    ,log10_pval := -log10(p)
    ]%>%inner_join(tfbs_cluster,by=c('elements'= 'Name'))%>%
    dplyr::select(-'numbtfbs_snps_notin_clust')%>%
    mutate('log10_numbsnps'=log10(numbtfbs_snps_in_clust))%>%
    mutate('log10_odds_ratio'=log10(odds_ratio+0.1))%>%
    dplyr::select(c('elements','ancestry','significance','elements','log10_pval','log10_numbsnps','log10_odds_ratio'))
  return(df)
}

denisova_clust_OR = reshape_OR(denisova_clust_OR,tfbs_cluster[[1]])
neanderthal_clust_OR = reshape_OR(neanderthal_clust_OR,tfbs_cluster[[2]])

denisova_highfreq_tfbs_OR = reshape_OR(denisova_highfreq_clust_OR,tfbs_cluster[[1]])
neanderthal_highfreq_tfbs_OR = reshape_OR(neanderthal_highfreq_clust_OR,tfbs_cluster[[2]])


## combine archaic enrichment and then plot it

asnps_clust_OR = rbind(denisova_clust_OR,neanderthal_clust_OR)
highfreq_asnps_clust_OR = rbind(denisova_highfreq_tfbs_OR,neanderthal_highfreq_tfbs_OR)

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per TF':= log10_numbsnps]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
  text=ifelse(df$significance!='',df$elements,'')
  
  ggplot(df,aes(x=log10_odds_ratio,y=log10_pval,label = text,col=log10_numbsnps))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
    geom_hline(yintercept=1.3, linetype="dashed", color = "black",size=0.2)+ ## threshold pval > 0.5 
    geom_text_repel(
      size = 5,
      color='black',
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines"),
      max.overlaps = getOption("ggrepel.max.overlaps", default = 20)
    )+
    gradient+
    xlab('\n log10 Odds ratio \n')+
    ylab('-Log10 (P)')+
    # xlim(-1,2)+
    facet_wrap(ancestry~.,ncol = 2)+
    theme(strip.text.x = element_text(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"),
          axis.line = element_blank()
          )
}


pdf(paste(plot_dir,'allsnps_tfbscluster_enrichment.pdf',sep=''),width = 10,height=6)
tf_enrich_plot(asnps_clust_OR)
dev.off()

pdf(paste(plot_dir,'highfreq_snps_tfbscluster_enrichment.pdf',sep=''),width = 10,height=6)
tf_enrich_plot(highfreq_asnps_clust_OR)
dev.off()



