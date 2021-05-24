## use this script to analyse the results of MotifbreakR
## in particular look at:
## 1. Ts:Tv ratio across each tisssue
## 2. Predicted impact on PMWs
## 3. Enrichment over TF clusters (Vierstra et al. 2020)
## 4. Jaccard score between sets of TFs affected by each ancestry

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

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'

tfbs_input_dir = './Motifbreak/output_files'
snps_input_dir = './Chromatin_states/SNPs_chromHMM_annotated'

columns_to_read = c(1,7:12,4,5,15)

immune_cells = c('BCells','TCells')

##-------------
## read files
##-------------
## cres snps
cres_snps = read_cres(snps_input_dir,columns_to_read)
names(cres_snps) = gsub("\\.txt$","",list.files(snps_input_dir,recursive = F,full.names = F))

## snps with tfbs predictions
hocomoco_tfbs = read_tfbs(tfbs_input_dir,'hocomoco')
lapply(hocomoco_tfbs,function(x)x[,c(..range_keys)] %>% unique() %>% nrow())

jaspar_tfbs = read_tfbs(tfbs_input_dir,'jaspar')
lapply(jaspar_tfbs,function(x)x[,c(..range_keys)] %>% unique() %>% nrow())

## combine hocomoco and jaspar results
combined_tfbs = purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
  lapply(
    function(x)x=x[
      ,alleleDiff:=round(alleleDiff,1)
    ][
      ,c(..range_keys,'REF','ALT','geneSymbol','providerName','alleleDiff')
      ] %>% unique()
)
combined_tfbs=lapply(combined_tfbs,function(x)x=x[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start,end,REF,ALT,geneSymbol)])

## combine tfbs with cell-type/state info
tfbs_cres = purrr::map2(cres_snps,combined_tfbs,inner_join,by=c(range_keys,'REF',"ALT"))
tfbs_cres = Map(mutate,tfbs_cres,'pop'=names(tfbs_cres))
tfbs_cres = lapply(tfbs_cres,function(x)x=x[MAF>=0.2])


## get some %
numb_tfbs_snps_cres = lapply(tfbs_cres,function(x)count_snps(x))
numb_cres_snps = lapply(cres_snps,function(x)count_snps(x[MAF>=0.2]))

prop_tfbs_cres = purrr::map2(numb_tfbs_snps_cres,numb_cres_snps,function(x,y)round(x/y*100,1))
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
        ][
          ,c(..range_keys,'REF','ALT','tv_ts','cell_type','cell_line')
          ]%>%unique()
  y=y[
      ,numb_ts:=sum(tv_ts =='ts'),by=.(cell_type)
      ][
        ,numb_tv:=sum(tv_ts =='tv'),by=.(cell_type)
        ][
          ,ts_tv_ratio:=round(numb_ts/numb_tv,2)
          ]%>%unique()
  y=y[
    ,c('ts_tv_ratio','cell_type','cell_line')
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
all_tfbs_cres_tstv = lapply(tfbs_cres,function(x)x=tvts_ratio(x)) 

## Make Ts:Tv ratio for all cells and report significance using one sample wilcoxon test
## from this remove the cell lines with less than 3 cell types
## test whether aSNPs Ts:Tv ratio is significantly different from that of naSNPs 
all_tfbs_cres_tstv =lapply(
  all_tfbs_cres_tstv,function(x)x=x[
    ,numb_celltypes:=.N,by=.(cell_line)
    ][
      numb_celltypes > 3
      ][
        ,numb_celltypes:=NULL
      ]
)
all_tfbs_cres_tstv = Map(mutate,all_tfbs_cres_tstv,pop=names(all_tfbs_cres_tstv))%>%rbindlist()

## function to create summary stat significance table
get_signif_table = function(x){
  df= copy(x)[
    ,p.signif:=ifelse(p_signif==' ','',p_signif)
    ][
      ,c('cell_line','p','p_adj','p_signif','group2')
      ]%>%unique()%>%setnames(old='group2',new='pop')
  
  ## add fake modern humans columns (otherwise u cant merge the results)
  fake_modern_stat =data.table(
    cell_line=df$cell_line,
    p = rep(0,nrow(df)),
    p_adj = rep(0,nrow(df)),
    p_signif = rep('',nrow(df)),
    pop = rep('modern_humans',nrow(df))
  )%>%unique()
  
  df =rbind(df,fake_modern_stat)%>%setorderv('cell_line',1)
  return(df)

}

stat_test_tstv = compare_means(
  ts_tv_ratio~pop,
  all_tfbs_cres_tstv,
  method = "wilcox.test",
  paired = FALSE,
  group.by = 'cell_line',
  ref.group = 'modern_humans'
) %>% as.data.table()%>%adjust_pvalues()%>%get_signif_table()

all_tfbs_cres_tstv = merge(all_tfbs_cres_tstv,stat_test_tstv,by=c('cell_line','pop'))

plot_tfbs = function(x,cells,y_scale,column,yintercept){
  df=copy(x)[
    ,column_to_plot:= column
    ][
          cell_line %in% c(cells)
          ]
  pvals=copy(df)[,c('pop','cell_line','p_signif')] %>% unique()
  modern_humans = copy(df)[pop %like% 'modern']%>%group_by(cell_line) %>%summarize(int = median(column_to_plot))%>%as.data.table()
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  
  ggplot(df,aes(x=pop,y=column_to_plot,fill=pop))+
    scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Modern  humans','Neanderthal'))+
    xlab(' ')+ylab(' ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_text(data=pvals, aes(x=pop, y=max(df$column_to_plot)+0.1, label=p_signif), col='black', size=7)+
    geom_hline(yintercept=yintercept, linetype="dashed", color = "darkgrey",size=0.2)+
    facet_wrap(cell_line~., ncol = 4,scales = y_scale)+
    geom_hline(data = modern_humans, aes(yintercept=int), linetype="dashed", color = "lightgrey",size=0.2)+
    theme(
      strip.text.x = element_text(),
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
      axis.line = element_line(color = "black", size = 0, linetype = "solid")
        )
}

pdf(paste(plot_dir,'tv_ts_ratio_all_cells.pdf',sep=''),width = 20,height = 15)
plot_tfbs(
  all_tfbs_cres_tstv,
  all_tfbs_cres_tstv$cell_line,
  'free_y',
  all_tfbs_cres_tstv$ts_tv_ratio,
  yintercept=2
  )
dev.off()

pdf(paste(plot_dir,'tv_ts_ratio_imunue_cells.pdf',sep=''),width = 8,height = 3)
plot_tfbs(
  all_tfbs_cres_tstv,
  immune_cells,
  'fixed',
  all_tfbs_cres_tstv$ts_tv_ratio,
  yintercept=2
  )
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
    ]%>%unique()
)%>%rbindlist()

stat_test_pwm_impact = compare_means(
  alleleDiff~pop,
  tfbs_cres_impact,
  method = "wilcox.test",
  paired = FALSE,
  group.by = 'cell_line',
  ref.group = 'modern_humans'
) %>% as.data.table()%>%adjust_pvalues()%>%get_signif_table()

tfbs_cres_impact = merge(tfbs_cres_impact,stat_test_pwm_impact,by=c('cell_line','pop'))

pdf(paste(plot_dir,'TFBS_effect_by_tissue.pdf',sep=''),width = 20,height = 20)
plot_tfbs(
  tfbs_cres_impact,
  tfbs_cres_impact$cell_line,
  'fixed',
  tfbs_cres_impact$alleleDiff,
  yintercept=0
  )
dev.off()

pdf(paste(plot_dir,'TFBS_effect_immune_cells.pdf',sep=''),width = 8,height = 3)
plot_tfbs(
  x = tfbs_cres_impact,
  cells = immune_cells,
  y_scale = 'fixed',
  tfbs_cres_impact$alleleDiff,
  yintercept=0
  )
dev.off()

##----------------------------------------------------------
## look enrichment across tf clusters (Vierstra et al 2020)
##----------------------------------------------------------
## first look enrichment for all TFBS-SNPs and then only for those with MIAF/DAF >= 0.2

motif=fread('../Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
]
# count snps in clusters
# NB: not all motifs get assigned because some might differ between pwm versions
snps_in_cluster = function(x){
  df = copy(tfbs_cres[[3]])%>%
    inner_join(motif,by=c('providerName'='motif_id'))%>%
      dplyr::select(c(all_of(range_keys),'Name','Cluster'))%>%unique()
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


## because some clusters might not be present in one of the 2 datasets and u are doing a full join
## some odds ratios might have the upper/lower ci = Inf
## set these to 0
denisova_clust_OR = calculate_odds_ratio(tfbs_cluster[[1]],tfbs_cluster[[2]],'Name',T)[,ancestry:='Denisova']%>%adjust_pvalues()
neanderthal_clust_OR = calculate_odds_ratio(tfbs_cluster[[3]],tfbs_cluster[[2]],'Name',T)[,ancestry:='Neanderthal']%>%adjust_pvalues()

## now prepare data.table for plotting
reshape_OR = function(x,y){
  df = copy(x)%>% mutate_if(is.numeric, function(y) ifelse(is.infinite(y), 1, y))%>%as.data.table()
  df = df[
    ,log10_padj := -log10(p_adj)
    ]%>%inner_join(y,by=c('elements'= 'Name'))%>%
    dplyr::select(-'numbtfbs_snps_notin_clust')%>%
    mutate('log10_numbsnps'=log10(numbtfbs_snps_in_clust))%>%
    mutate('log10_odds_ratio'=log10(odds_ratio+0.1))%>%
    dplyr::select(c('elements','ancestry','significance','elements','p','p_adj','log10_padj','log10_numbsnps','log10_odds_ratio'))
  return(df)
}

denisova_clust_OR = reshape_OR(denisova_clust_OR,tfbs_cluster[[1]])
neanderthal_clust_OR = reshape_OR(neanderthal_clust_OR,tfbs_cluster[[3]])

## combine archaic enrichment and then plot it

asnps_clust_OR = rbind(denisova_clust_OR,neanderthal_clust_OR)

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per Cluster':= log10_numbsnps]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per Cluster`),option="inferno",discrete = F)
  text=ifelse(df$significance!='',df$elements,'')
  
  ggplot(df,aes(x=log10_odds_ratio,y=log10_padj,label = text,col=log10_numbsnps))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
    geom_hline(yintercept=1.3, linetype="dashed", color = "black",size=0.2)+ ## threshold pval > 0.05 
    geom_text_repel(
      size = 5,
      color='black',
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines"),
      max.overlaps = getOption("ggrepel.max.overlaps", default = 20)
    )+
    gradient+
    xlab('\n log10 Odds ratio \n')+
    ylab('-Log10 (P adj)')+
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

pdf(paste(plot_dir,'allsnps_tfbscluster_enrichment.pdf',sep=''),width = 15,height=10)
tf_enrich_plot(asnps_clust_OR)
dev.off()


## calculate jaccard index
jaccard_tfbs_cres =copy(tfbs_cres)%>%lapply(
  function(x)x=split(x,by='cell_type')%>%
  lapply(
    function(y)y=y[,geneSymbol]%>%unique()
    )
)


jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}


jaccard_score_denisova = purrr::map2(jaccard_tfbs_cres[[1]],jaccard_tfbs_cres[[2]],function(x,y) jaccard(x,y)%>%data.table(score=.))
jaccard_score_denisova = Map(mutate,jaccard_score_denisova,cell_type=names(jaccard_score_denisova))%>%rbindlist()%>%group_celltypes()
jaccard_score_denisova =jaccard_score_denisova[,pop:='denisova'][,mean_score:=mean(score)]


jaccard_score_neanderthal = purrr::map2(jaccard_tfbs_cres[[3]],jaccard_tfbs_cres[[2]],function(x,y) jaccard(x,y)%>%data.table(score=.))
jaccard_score_neanderthal = Map(mutate,jaccard_score_neanderthal,cell_type=names(jaccard_score_neanderthal))%>%rbindlist()%>%group_celltypes()
jaccard_score_neanderthal =jaccard_score_neanderthal[,pop:='neanderthal'][,mean_score:=mean(score)]

jaccard_score_archaics =rbind(jaccard_score_denisova,jaccard_score_neanderthal)
jaccard_score_archaics = jaccard_score_archaics[,numbcells:=.N,by=.(cell_line,pop)][numbcells>3]

stat_test_jaccard = compare_means(
  score~pop,
  jaccard_score_archaics,
  method = "wilcox.test",
  paired = FALSE,
  group.by = 'cell_line',
  ref.group='denisova'
) %>% as.data.table()%>%adjust_pvalues()


jaccard_plot = ggplot(jaccard_score_archaics,aes(x=pop,y=score,fill=pop))+
    scale_fill_manual(name= " ",values = my_palette_pop[-2],labels = c('Denisova','Neanderthal'))+
    xlab(' ')+ylab(' Jaccard similarity score ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    facet_wrap(cell_line~., ncol = 4,scales = 'fixed')+
    theme(
      strip.text.x = element_text(),
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
      axis.line = element_line(color = "black", size = 0, linetype = "solid")
)

pdf(paste(plot_dir,'jaccard_similarity_allcells.pdf',sep=''),width = 10,height=10)
jaccard_plot
dev.off()


## collect and write all pval tables and make a single excel file

table_pwm_impact = copy(stat_test_pwm_impact)[pop != 'modern_humans']%>%setnames(old='cell_line',new='tissue')

table_tstv = copy(stat_test_tstv)[pop != 'modern_humans']%>%setnames(old='cell_line',new='tissue')

table_jaccard = copy(stat_test_jaccard)[
  ,comparison:=paste(group1,group2,sep=' vs ')
  ][
    ,c('cell_line','p','p_adj','p_signif','comparison')
]%>%setnames(old='cell_line',new='tissue')

spreadsheet = list(table_pwm_impact,table_tstv,table_jaccard)
names(spreadsheet) = c('Impact on PWMs','TsTv','Jaccard similarity')

write.xlsx(spreadsheet,paste(table_dir,'Supp_Table_TFBS_pvals.xlsx',sep=''),append=T)
