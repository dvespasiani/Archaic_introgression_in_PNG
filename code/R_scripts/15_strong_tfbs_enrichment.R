## look at proportion STRONG effect tfbs SNPs across cells
library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)

options(width=150)
numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'

motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'

merging_keys=c('seqnames','start','end')
immune_cells=c('TCells')

cell_line_levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
                  'iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus')

allele_frequency=0.3

##-----------------------------------------
##         Recurrent function
##-----------------------------------------

assign_names=function(x){
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}



##-----------------------------------------
##         Read files
##-----------------------------------------
#1)
read_tfbs_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,end:=start+1
        ])
  x=assign_names(x)
  return(x)
}


tfbs=read_tfbs_snps(paste(motif_input_dir,'new_set/combined/',sep=''))
lapply(tfbs,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

#2)
read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq'))[
        chrom_state%in%c('2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
        ][
          ,chrom_state:=NULL
          ])
  x=assign_names(x)
  return(x)
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') 

## combine tfbs with cell-type/state info
tfbs_cres=purrr::map2(tfbs,states,inner_join,by=merging_keys)
tfbs_cres=Map(mutate,tfbs_cres,'pop'=names(tfbs_cres))%>% 
  lapply(function(x)x=as.data.table(x)[
    ,c('seqnames','start','end','REF','ALT','cell_type','cell_line','all_freq','effect','alleleDiff','pop')
    ] %>% unique())

lapply(tfbs_cres,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

asnps_enrichment=function(x,pwm_effect){

  calculate_enrichment=function(x){
    df=copy(x)
    df=lapply(df,function(y)y=y[
        ,numb_tfbs_snps_cell_type:=.N, by =.(cell_type,cell_line)
        ][
          ,c('numb_tfbs_snps_cell_type','cell_type','cell_line','pop')] %>% unique())

    deni_enrichment=inner_join(df[[1]],df[[3]],by=c('cell_type','cell_line')) %>% as.data.table()
    deni_enrichment=deni_enrichment[
      ,asnps_nasnps_ratio:=numb_tfbs_snps_cell_type.x/numb_tfbs_snps_cell_type.y
      ][
        ,mean_asnps_nasnps_ratio:=mean(asnps_nasnps_ratio)
        ][
        ,log2_enrichm:=log2(((numb_tfbs_snps_cell_type.x/numb_tfbs_snps_cell_type.y)/(mean_asnps_nasnps_ratio)))
        ]
    nean_enrichment=inner_join(df[[2]],df[[3]],by=c('cell_type','cell_line')) %>% as.data.table()
    nean_enrichment=nean_enrichment[
      ,asnps_nasnps_ratio:=numb_tfbs_snps_cell_type.x/numb_tfbs_snps_cell_type.y
      ][
        ,mean_asnps_nasnps_ratio:=mean(asnps_nasnps_ratio)
        ][
          ,log2_enrichm:=log2(((numb_tfbs_snps_cell_type.x/numb_tfbs_snps_cell_type.y)/(mean_asnps_nasnps_ratio)))
          ]
    asnps_tfbs_enrichment=list(deni_enrichment,nean_enrichment)
    return(asnps_tfbs_enrichment)
  }

  df=copy(x)
  df=lapply(df,function(y)y=y[,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][effect%in%pwm_effect][,c('effect'):=NULL] %>% unique())
  df=calculate_enrichment(df)
  df=lapply(df,function(a)a=a[,c('pop.x','cell_type','cell_line','log2_enrichm')] %>% 
              setnames(old = 'pop.x',new='pop')%>% unique()) %>% rbindlist()
  return(df)

}

tfbs_asnps_enrichment=asnps_enrichment(tfbs_cres,'strong')

# lapply(tfbs_cres,function(x)x[effect%in%'strong'][,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][,c('seqnames','start','end')] %>% unique() %>% nrow())
# lapply(tfbs_cres,function(x)x[effect%in%'strong'][cell_line%in%'TCells'][,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][,c('seqnames','start','end')] %>% unique() %>% nrow())
# lapply(tfbs_cres,function(x)x[effect%in%'strong'][cell_line%in%'BCells'][,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][,c('seqnames','start','end')] %>% unique() %>% nrow())
lapply(tfbs_cres,function(x)x[effect%in%'strong'][cell_line%in%c('BCells','TCells')][,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][,c('seqnames','start','end')] %>% unique() %>% nrow())

## calculate significance via one-tailed t test
stat_test_enrich=copy(tfbs_asnps_enrichment)[
  !cell_line%in%c('Adipose','IMR90_fetal_lung_fibroblast','Myosatellite','Thymus') ## remove these cell line as they lack replicates
] %>% unique()
stat_test_enrich=stat_test_enrich %>% split(as.factor(stat_test_enrich$cell_line))

one_sample_ttest=function(x){
  df=copy(x)
  deni=copy(df)[pop%in%'denisova']
  nean=copy(df)[pop%in%'neandertal']
  
  perform_test=function(x){
    df=copy(x)
    df=df[
      ,p:=t.test(log2_enrichm,mu=0)$p.val
        ][
          ,statistic:=t.test(log2_enrichm,mu=0)$stat
          ][
            ,df:=t.test(log2_enrichm,mu=0)$parameter
            ]
      return(df)
    
  }
  
  deni=perform_test(deni)
  nean=perform_test(nean) 
  combined=rbind(deni,nean) 
  return(combined)
}

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

stat_test_enrich=lapply(stat_test_enrich,function(x)one_sample_ttest(x)) %>% rbindlist() %>% adjust_pvalues()

tfbs_asnps_enrichment$cell_line=factor(tfbs_asnps_enrichment$cell_line,
levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
         'iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus'))


### write supp table
enrich_pval_tables=function(x,table){
  df_pval=copy(x)[,c('pop','cell_line','statistic','df','p','p.adj','p.signif')] %>% unique()
  df_nonpval=copy(x)[,c('statistic','p','p.adj','p.signif','df'):=NULL] %>% unique()
  
  dfs=list(df_nonpval,df_pval)
  
  df_names=c(table,'pvalues') 
  names(dfs)=df_names
  for (i in seq_along(dfs)){
    assign(df_names[i],dfs[[i]],.GlobalEnv)}
  return(dfs) 
}


tfbs_tables=enrich_pval_tables(stat_test_enrich,'fraction_snps_per_cell_type')
write.xlsx(tfbs_tables,paste(table_dir,'Supp_Table_7_TFBS_pvalues.xlsx',sep=''))

tfbs_asnps_enrichment_plot=function(df,cells,y_scale){
  
  df=df[cell_line%in%cells]
  df$cell_line=factor(df$cell_line,levels=cell_line_levels)
  stat_test_prop=stat_test_enrich[cell_line%in%cells][,c('pop','cell_line','p.signif')] %>% unique()
  stat_test_prop$cell_line=factor(stat_test_prop$cell_line,levels=cell_line_levels)
  
  
  my_palette_pop=c('#C99E10', # denisova
    '#9B4F0F' #neandertal
    )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal'))
  
  ggplot(df,aes(x=pop,y=log2_enrichm,fill=pop))+
    geom_violin(trim=F,scale = "width")+colScale+
    xlab(' ')+ylab('\n Log2 aSNPs enrichment \n ')+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_jitter(position=position_jitter(0.1),size=0.2)+
    geom_hline(yintercept = 0, linetype="dashed", color = "black",size=0.2)+
    geom_text(data=stat_test_prop, aes(x=pop, y=max(df$log2_enrichm)+0.1, label=p.signif), col='black', size=7)+
    facet_wrap(cell_line~., ncol = 4,scales =y_scale)+
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

pdf(paste(plot_dir,'tfbs_asnps_enrichment_across_all_cells.pdf',sep=''),width = 20,height = 20)
tfbs_asnps_enrichment_plot(tfbs_asnps_enrichment,tfbs_asnps_enrichment$cell_line,'free_y')
dev.off()

pdf(paste(plot_dir,'tfbs_asnps_enrichment_immune_cells.pdf'),width = 8,height = 4)
tfbs_asnps_enrichment_plot(tfbs_asnps_enrichment,immune_cells,'fixed')
dev.off()

##-----------------------------------------
##         Compute Ts:Tv ratio 
##-----------------------------------------
tvts_ratio=function(x,pwm_effect){
  x=copy(x)[,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][effect%in%pwm_effect]
  y=copy(x)[
    ,dna_ref_base:=ifelse(REF=='A' |REF=='G','purine','pyrimidine')
    ][
      ,dna_alt_base:=ifelse(ALT=='A' |ALT=='G','purine','pyrimidine')
      ][
        ,tv_ts:=ifelse(dna_ref_base=='purine' & dna_alt_base=='pyrimidine','tv',
                       ifelse(dna_ref_base=='pyrimidine' & dna_alt_base=='purine','tv','ts'))
        ]
  y=unique(y[
    ,c('seqnames','start','end','REF','ALT','tv_ts','pop','cell_type','cell_line')
    ][
      ,numb_ts:=sum(tv_ts =='ts'),by=.(cell_type)
      ][
        ,numb_tv:=sum(tv_ts =='tv'),by=.(cell_type)
        ][
          ,ts_tv_ratio:=numb_ts/numb_tv
          ])
  y=y[
    ,c('tv_ts','ts_tv_ratio','pop','cell_type','cell_line')
    ]%>% unique() 
  return(y)
}

cells_tstv_ratio=lapply(tfbs_cres,function(x)x=tvts_ratio(x,'strong')) %>% rbindlist()

## Make Ts:Tv ratio for all cells and report significance using one sample wilcoxon test
stat_test=copy(cells_tstv_ratio)
stat_test=stat_test %>% split(as.factor(stat_test$pop)) %>%
  lapply(function(x)x=x %>% split(as.factor(x$cell_line)) %>%
           lapply(function(y)
             y=y %>% 
               mutate('p'=wilcox.test(y$ts_tv_ratio,mu=2,exact = F,paired = F)$p.val)%>%
               mutate('statistic'=wilcox.test(y$ts_tv_ratio,mu=2,exact = F,paired = F)$stat)%>% 
               as.data.table()
           )%>%rbindlist()
  ) %>%rbindlist() %>% adjust_pvalues()

## make two Supp tables (fraction snps and pvals)
enrich_pval_tables=function(x,table){
  df_pval=copy(x)[,c('pop','cell_line','statistic','p','p.adj','p.signif')] %>% unique()
  df_pval=df_pval[
    ,p:=ifelse(is.nan(p),0,p)
  ][
    ,p.adj:=ifelse(is.nan(p.adj),0,p.adj)
    ][
      ,p.signif:=ifelse(is.na(p.signif),' ',p.signif)
      ]
  df_nonpval=copy(x)[,c('statistic','p','p.adj','p.signif'):=NULL] %>% unique()
  
  dfs=list(df_nonpval,df_pval)
  
  df_names=c(table,'pvalues') 
  names(dfs)=df_names
  for (i in seq_along(dfs)){
    assign(df_names[i],dfs[[i]],.GlobalEnv)}
  return(dfs) 
}


ts_tv_tables=enrich_pval_tables(stat_test,'Ts_Tv_ratio')
write.xlsx(ts_tv_tables,paste(table_dir,'Supp_Table_8_TsTv_pvalues.xlsx',sep=''))

## average png ts:tv ratio
png_tstv=copy(ts_tv_tables[[1]])[pop=='png'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique() ## 1.9

denisova_tstv_immune=copy(ts_tv_tables[[1]])[pop=='denisova'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique()
nean_tstv_immune=copy(ts_tv_tables[[1]])[pop=='neandertal'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique() 
png_tstv_immune=copy(ts_tv_tables[[1]])[pop=='png'][cell_line%in%'TCells'][,mean:=mean(ts_tv_ratio)][,mean] %>% unique() 


cells_tstv_ratio$cell_line=factor(cells_tstv_ratio$cell_line,
levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
         'iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus'))

tstv_ratio_plot=function(df,cells,y_scale){
  df=df[cell_line%in%cells]
  df$cell_line=factor(df$cell_line,levels=cell_line_levels)
  
  stat_test=stat_test[cell_line%in%cells][,c('pop','cell_line','p.signif')] %>% unique()
  stat_test$cell_line=factor(stat_test$cell_line,levels=cell_line_levels)
  
  
  png=copy(df)[pop%in%'png']
  png=png%>%group_by(cell_line) %>%summarize(int = median(ts_tv_ratio))
  
  
  my_palette_pop=c(#'#8D230F', #ambiguous
    '#C99E10', # denisova
    '#9B4F0F', #neandertal
    '#1E656D'#png
  )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
  
  ggplot(df,aes(x=pop,y=ts_tv_ratio,fill=pop))+
    colScale+
    xlab(' ')+ylab('\n Ts:Tv ratio \n ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_jitter(position=position_jitter(0.1),size=0.2)+
    geom_text(data=stat_test, aes(x=pop, y=3.05, label=p.signif), col='black', size=7)+
    geom_hline(yintercept=2, linetype="dashed", color = "darkgrey",size=0.2)+
    facet_wrap(cell_line~., ncol = 4,scales = y_scale)+
    geom_hline(data = png, aes(yintercept=int), linetype="dashed", color = "lightgrey",size=0.2)+
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

pdf(paste(plot_dir,'tv_ts_ratio_all_cells.pdf'),width = 20,height = 20)
tstv_ratio_plot(cells_tstv_ratio,cells_tstv_ratio$cell_line,'free_y')
dev.off()

pdf(paste(plot_dir,'tv_ts_ratio_imunue_cells.pdf'),width = 8,height = 4)
tstv_ratio_plot(cells_tstv_ratio,immune_cells,'fixed')
dev.off()

# ##----------------------------
# ##     ATAC-seq peaks
# ##----------------------------
# 
# read_atacseq=function(path){
#   state=sub('.*\\/', '',path)
#   state=sub("\\_.*","",state)
#   x=fread(path,sep='\t',header = T)[
#     ,c('seqnames','start','end'):= tstrsplit(`peak_id` , "_", fixed=TRUE)
#     ][
#       ,c('seqnames','start','end')
#       ][
#         ,start:=as.numeric(start)
#         ][
#           ,end:=as.numeric(end)
#           ][
#             ,cell_state:=state
#             ] %>% setorderv(c('seqnames','start'),1)%>% unique()
#   return(x)
# 
# }
# 
# rested_atacseq=read_atacseq('../Annotation_and_other_files/Calderon_atac_seq_peaks/rested_cells.gz')
# stimulated_atacseq=read_atacseq('../Annotation_and_other_files/Calderon_atac_seq_peaks/stimulated_cells.gz')
# immune_atacseq=rbind(rested_atacseq,stimulated_atacseq)
# 
# # ## Immune strong tfbs snps in atac peaks 
# selected_tfbs_atac_peaks=function(tfbs,cells){
#   tfbs=lapply(tfbs,function(x)x=x[cell_line%in%cells])
#   setkeyv(immune_atacseq,merging_keys)
#   lapply(tfbs,function(x)setkeyv(x,merging_keys))
# 
#   selected_snps_in_ataseq=lapply(tfbs,function(x)
#     x=foverlaps(x,immune_atacseq,type = 'within')[
#       ,c('seqnames','i.start','i.end','cell_state','pop','all_freq','effect')
#       ][
#         ,start:=i.start
#         ][
#           ,end:=i.end
#           ][
#             ,c('i.start','i.end'):=NULL
#             ] %>% na.omit()%>% unique()
#   )
#   return(selected_snps_in_ataseq)
# }
# immune_tfbs_atac=selected_tfbs_atac_peaks(new_tfbs_cres,immune_cells)
# # 
# # ## fraction in rested/stimulated
# stimulated_vs_rested=function(x,pwm_effect){
#   x=copy(x)
#   x=x %>% lapply(function(y)y=y[,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][
#     effect%in%pwm_effect])
#   x=x%>% lapply(function(y)
#     y=y[
#       ,'length':=.N,by=.(seqnames,start,end)
#       ][
#         ,'cell_condition':=ifelse(length==2,'both',`cell_state`)
#         ][
#           ,c('seqnames','start','end','cell_condition','pop')
#           ] %>% unique()
#   )
#   x=lapply(x,function(y)
#     y=y[
#       ,number_snps:=.N
#       ][
#         ,number_snps_per_condition:=.N,by=.(cell_condition)
#         ][
#           ,fraction:=number_snps_per_condition/number_snps
#           ][
#             ,c('number_snps','number_snps_per_condition','cell_condition','fraction','pop')
#             ] %>% unique()) %>% rbindlist()
#   x=x[
#     ,fraction:=round(fraction,2)
#     ]
# }
# 
# immune_stimulated_rested=stimulated_vs_rested(immune_tfbs_atac,'strong')
# ## enrichment in these elements
# open_chrom_regions=copy(immune_atacseq) %>% unique()
# open_chrom_regions=open_chrom_regions[
#   ,'length':=.N,by=.(seqnames,start,end)
#   ][
#     ,'cell_condition':=ifelse(length==2,'both',`cell_state`)
#     ][,cell_state:=NULL][
#       ,number_snps:=.N
#       ][
#         ,number_snps_per_condition:=.N,by=.(cell_condition)
#         ][
#           ,c('cell_condition','number_snps_per_condition','number_snps')
#           ][
#             ,fraction:=number_snps_per_condition/number_snps
#             ][
#               ,pop:='Calderon et al.'
#               ] %>% unique()
# 
# immune_stimulated_rested=rbind(immune_stimulated_rested,open_chrom_regions)
# 
# test_enrichment=function(file,condition,selected_pop){
#   df=copy(file)[pop%in% c('denisova','neandertal') & cell_condition%in%condition][
#     ,fraction:=NULL
#     ][
#       ,totsnps_percond:=sum(number_snps_per_condition)
#       ][
#         ,tot_other_asnps:=ifelse(selected_pop=='denisova',number_snps[[2]],number_snps[[1]])
#         ][
#           pop%in%selected_pop
#           ][
#             ,test:=phyper(number_snps_per_condition-1,number_snps,tot_other_asnps,totsnps_percond,lower.tail = F)
#             ]
# }
# 
# nean_enrichment=test_enrichment(immune_stimulated_rested,'rested','neandertal')
# 
# deni_enrichment=test_enrichment(immune_stimulated_rested,'stimulated','denisova')
# 
# 
# immune_stimulated_rested$pop=factor(immune_stimulated_rested$pop,levels=c('denisova','neandertal','png','Calderon et al.'))
# 
# ## plot 
# my_palette_cells=c( 'azure3', # both
#                     'ivory2',#rested
#                     'lightgoldenrod3'# stimulated
# )
# 
# names(my_palette_cells)=levels(as.factor(immune_stimulated_rested$cell_condition))
# colScale_cells=scale_fill_manual(name= " ",values = my_palette_cells,labels = c('Both','Rested','Stimulated'))
# 
# cell_condition_plot=ggplot(immune_stimulated_rested,aes(x=pop,y=fraction,fill=cell_condition))+
#   geom_bar(stat='identity')+
#   colScale_cells+
#   scale_x_discrete(labels=c("denisova" = "Denisova","neandertal" = "Neandertal",'png'='PNG','Calderon et al.'='Calderon et al.'))+
#   ylab('\n Fraction TFBS-SNPs in open chromatin regions  \n')+
#   xlab(' ')+
#   theme(panel.background =element_rect(fill = 'white',size = 1, colour = 'black'),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         legend.text = element_text(),
#         legend.title = element_text(),
#         axis.text.x = element_text(),
#         axis.text.y = element_text(),
#         axis.title.y = element_text(hjust=0.5),
#         axis.text=element_text(),
#         axis.line = element_line(color = "black", size = 0, linetype = "solid"))
# 
# 
# pdf('/home/dvespasiani/Archaic_introgression_in_PNG/tfbs_plots/atac_immune_tfbs_snps.pdf',width = 5,height =4)
# cell_condition_plot
# dev.off()

