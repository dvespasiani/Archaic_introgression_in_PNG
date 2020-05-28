## look at proportion STRONG effect tfbs SNPs across cells
library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)

setDTthreads(10)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')
plot_output_dir='./home/dvespasiani/tfbs_plots/'
merging_keys=c('seqnames','start','end')
immune_cells=c('TCells','BCells')


gtex_expr=fread('./GTEx/GTEx_Analysis_v8_gene_median_tmp/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',sep='\t',header = T)
gene_ids=gtex_expr[,c(1:2)]


tfbs=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            select =c('seqnames','start','end','REF','ALT','geneSymbol','effect')) %>% unique())
      x=lapply(x,function(y)
        y=y[
              ,duplicated_effect:=.N,by=.(seqnames,start,end,geneSymbol)
              ]%>% group_by(seqnames,start,end,geneSymbol) %>%
        filter(duplicated_effect==max(duplicated_effect)) %>% as.data.table()
        # [
        # ,delta_pct:=abs(pctRef-pctAlt)
        # ] %>% group_by(seqnames,start,end,geneSymbol) %>%
        # filter(delta_pct==max(delta_pct)) %>% as.data.table()
    )
  x=x[c(2:4)] # remove ambiguous like this for now
  x=lapply(x,function(y)y=semi_join(y,gene_ids,by=c('geneSymbol'='Description')) %>%
             as.data.table()) ## remove TFs absent from gtex (only D/OBOXs)
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  }

snps_tfbs=tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/combined/')


#LOOK AT PROPORTION OF STRONG TFBS ON ≠ CELLS FOR THE 3 PEOPLE
read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq'))[
              chrom_state%in%c('2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
              ][
                ,chrom_state:=NULL
                ] %>% unique()
    )
  # x=x[c(2:4)] # remove ambiguous like this for now
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') # remove the temp directory once the old script is optimized

snps_tfbs_cres=purrr::map2(snps_tfbs,states,inner_join,by=merging_keys)
snps_tfbs_cres=lapply(snps_tfbs_cres,function(x)setDT(x))

strong_highfreq=function(x){
  x=copy(x) %>% 
    lapply(function(y)
      y=y[all_freq>=0.2 & effect%in%'strong'])
  x=lapply(x,function(y)y=y[
    ,effect_strenght:=ifelse(effect=='weak',as.numeric(0),as.numeric(1))
    ] %>% group_by(seqnames,start,end) %>%
    filter(effect_strenght==max(effect_strenght)) %>% as.data.table() %>% unique()
  )
  x=lapply(x,function(y)
    y=y[
      ,number_snps_per_cell:=uniqueN(c(seqnames,start,end)),by=(cell_type)
      ][
        ,number_snps:=uniqueN(c(seqnames,start,end))
        
        ][
          ,fraction_snps_per_cell_per_effect:=number_snps_per_cell/number_snps
          ])
  # x=assign_names(x)
  x=Map(mutate,x,'pop'=names(x))
  x=lapply(x,function(y)setDT(y))
}


snps_tfbs_highfreq_strong=strong_highfreq(snps_tfbs_cres) 
# lapply(snps_tfbs_highfreq_strong,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
## 615 Amb
## 1802 Deni
## 1630 Nean
## 9936 PNG

selected_snps=copy(snps_tfbs_highfreq_strong)%>% rbindlist()
selected_snps=selected_snps[
  ,c('pop','fraction_snps_per_cell_per_effect','cell_type','cell_line')
  ] %>% unique()



stat_test_prop=copy(selected_snps)
stat_test_prop=stat_test_prop %>% split(as.factor(stat_test_prop$cell_line))

one_tailed_wilcoxon=function(x){
  df=copy(x)
  deni=copy(df)[pop%in%'denisova']
  nean=copy(df)[pop%in%'neandertal']
  png=copy(df)[pop%in%'png']
  
perform_test=function(x,y){
    x=x[,p:=wilcox.test(x$fraction_snps_per_cell_per_effect,y$fraction_snps_per_cell_per_effect,alternative = 'g',exact = F)$p.val
        ][
      ,statistic:=wilcox.test(x$fraction_snps_per_cell_per_effect,y$fraction_snps_per_cell_per_effect,alternative = 'g',exact = F)$stat
      ][
        ,p.adj:=p.adjust(p,method = 'bonferroni')
        ][
          ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                             ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                    ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                           ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
          ]
}

  deni=perform_test(deni,png)
  nean=perform_test(nean,png) 
  combined=rbind(deni,nean) 
  return(combined)
}

stat_test_prop=lapply(stat_test_prop,function(x)one_tailed_wilcoxon(x)) %>% rbindlist()


## make two Supp tables (fraction snps and pvals)
enrich_pval_tables=function(x,table){
  df_pval=copy(x)[,c('pop','cell_line','statistic','p','p.adj','p.signif')] %>% unique()
  df_nonpval=copy(x)[,c('statistic','p','p.adj','p.signif'):=NULL] %>% unique()
  
  dfs=list(df_nonpval,df_pval)

  df_names=c(table,'pvalues') 
    names(dfs)=df_names
    for (i in seq_along(dfs)){
      assign(df_names[i],dfs[[i]],.GlobalEnv)}
    return(dfs) 
  }



# tfbs_tables=enrich_pval_tables(stat_test_prop,'fraction_snps_per_cell_type')
# write.xlsx(tfbs_tables,'~/pvalue_tables/Supp_Table_8_TFBS_pvalues.xlsx')


selected_snps$cell_line=factor(selected_snps$cell_line,levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
                                                     'iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus'))


## 
selected_snps_plot=function(df,cells){
  
  df=df[cell_line%in%cells]
 
   stat_test_prop=stat_test_prop[cell_line%in%cells][,c('pop','cell_line','p.signif')] %>% unique()
  
      my_palette_pop=c(#'#8D230F', #ambiguous
                     '#C99E10', # denisova
                     '#9B4F0F', #neandertal
                     '#1E656D' #png
    )
    
names(my_palette_pop)= levels(as.factor(df$pop))
colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
  
  ggplot(df,aes(x=pop,y=fraction_snps_per_cell_per_effect,fill=pop))+
  geom_violin(trim=F,scale = "width")+colScale+
  xlab(' ')+ylab('\n Fraction TFBS-SNPs per epigenome \n ')+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_jitter(position=position_jitter(0.1),size=0.2)+
    geom_text(data=stat_test_prop, aes(x=pop, y=max(df$fraction_snps_per_cell_per_effect)+0.047, label=p.signif), col='black', size=7)+
   #   stat_compare_means(method = "wilcox.test",label = "p.signif",
  #                      ref.group = "png",hide.ns =T,
  #                      label.y =max(df$fraction_snps_per_cell_per_effect+0.05),
  #                      size=10)+
  facet_wrap(cell_line~., ncol = 4)+
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

pdf('/home/dvespasiani/tfbs_plots/prop_tfbs_snps_across_all_cells.pdf',width = 20,height = 20)
selected_snps_plot(selected_snps,selected_snps$cell_line)
dev.off()

pdf('/home/dvespasiani/tfbs_plots/prop_tfbs_snps_immune_cells.pdf',width = 7,height = 3)
selected_snps_plot(selected_snps,immune_cells)
dev.off()

## For these SNPs compute Ts:Tv ratio 
tvts_ratio=function(x){
  x=setDT(x)
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

cells_tstv_ratio=lapply(snps_tfbs_highfreq_strong,function(x)
 x=tvts_ratio(x)) %>% rbindlist()

## Make Ts:Tv ratio for all cells and report significance using one sample wilcoxon test
stat_test=copy(cells_tstv_ratio)
stat_test=stat_test %>% split(as.factor(stat_test$pop)) %>%
  lapply(function(x)x=x %>% split(as.factor(x$cell_line)) %>%
           lapply(function(y)
             y=y %>% 
               mutate('p'=wilcox.test(y$ts_tv_ratio,alternative = 'l',mu=2,exact = F,paired = F)$p.val)%>%
               mutate('statistic'=wilcox.test(y$ts_tv_ratio,alternative = 'l',mu=2,exact = F,paired = F)$stat)%>% 
               mutate('p.adj'=p.adjust(p,method = 'bonferroni'))%>%
               mutate('p.signif'=ifelse(`p.adj`<=0.0001,'****',
                                                 ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                                        ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                                               ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))) %>% 
               as.data.table()
           )%>%rbindlist()
  ) %>%rbindlist()

stat_test=stat_test[
  ,c('pop','cell_type','cell_line','ts_tv_ratio','statistic','p','p.adj','p.signif')
  ] %>% unique()

## make two Supp tables (fraction snps and pvals)

enrich_pval_tables=function(x,table){
  df_pval=copy(x)[,c('pop','cell_line','statistic','p','p.adj','p.signif')] %>% unique()
  df_nonpval=copy(x)[,c('statistic','p','p.adj','p.signif'):=NULL] %>% unique()
  
  dfs=list(df_nonpval,df_pval)
  
  df_names=c(table,'pvalues') 
  names(dfs)=df_names
  for (i in seq_along(dfs)){
    assign(df_names[i],dfs[[i]],.GlobalEnv)}
  return(dfs) 
}



# ts_tv_tables=enrich_pval_tables(stat_test,'Ts_Tv_ratio')
# write.xlsx(ts_tv_tables,'~/pvalue_tables/Supp_Table_9_TsTv_pvalues.xlsx')


tstv_ratio_plot=function(df,cells){
  df=df[cell_line%in%cells]
  stat_test=stat_test[cell_line%in%cells][,c('pop','cell_line','p.signif')] %>% unique()
  
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
    facet_wrap(cell_line~., ncol = 4)+
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

pdf('/home/dvespasiani/tfbs_plots/tv_ts_ratio_all_cells.pdf',width = 20,height = 20)
tstv_ratio_plot(cells_tstv_ratio,cells_tstv_ratio$cell_line)
dev.off()

pdf('/home/dvespasiani/tfbs_plots/tv_ts_ratio_imunue_cells.pdf',width = 7,height = 3)
tstv_ratio_plot(cells_tstv_ratio,immune_cells)
dev.off()

## ATAC-seq peaks
## calderon atac seq see if tfbs snps are more in stimulated/rested t/b immune cells 
read_atacseq=function(path){
 state=sub('.*\\/', '',path)
  state=sub("\\_.*","",state)
  x=fread(path,sep='\t',header = T)[
    ,c('seqnames','start','end'):= tstrsplit(`peak_id` , "_", fixed=TRUE)
    ][
      ,c('seqnames','start','end','contrast')
      ][
        ,start:=as.numeric(start)
        ][
          ,end:=as.numeric(end)
          ][
            ,cell_state:=state
          ]
  return(x)
  
}

rested_atacseq=read_atacseq('../Calderon_atac_seq_peaks/rested_cells')
stimulated_atacseq=read_atacseq('../Calderon_atac_seq_peaks/stimulated_cells')
immune_atacseq=rbind(rested_atacseq,stimulated_atacseq)

## Immune strong tfbs snps in atac peaks 
selected_tfbs_atac_peaks=function(tfbs,cells){
  tfbs=lapply(tfbs,function(x)x=x[cell_line%in%cells])
  setkeyv(immune_atacseq,merging_keys)
  lapply(tfbs,function(x)setkeyv(x,merging_keys))
  
 selected_snps_in_ataseq=lapply(tfbs,function(x)
    x=foverlaps(x,immune_atacseq,type = 'within')[
      ,c('seqnames','i.start','i.end','cell_state','pop')
      ][
        ,start:=i.start
        ][
          ,end:=i.end
          ][
            ,c('i.start','i.end'):=NULL
            ] %>% na.omit()%>% unique()
  )
  return(selected_snps_in_ataseq)
}
immune_tfbs_atac=selected_tfbs_atac_peaks(snps_tfbs_highfreq_strong,immune_cells)

# lapply(immune_tfbs_atac,function(x) x=x[,c('seqnames','start','end')] %>% unique() %>% nrow())
## 1446 Denisova SNPs (141 strong)
## 749 Neandertal (107 strong)
## 11024 PNG (522 strong)

## fraction in rested/stimulated
## if variant is present in both condition select the stimulated one 
stimulated_vs_rested=function(x){
  x=copy(x) %>% lapply(function(y)
  y=y[
    ,'length':=.N,by=.(seqnames,start,end)
    ][
    ,'cell_condition':=ifelse(length==2,'stimulated',`cell_state`)
    ][
      ,c('seqnames','start','end','cell_condition','pop')
    ] %>% unique()
)
x=lapply(x,function(y)
  y=y[
    ,number_snps:=.N
  ][
    ,number_snps_per_condition:=.N,by=.(cell_condition)
  ][
    ,fraction:=number_snps_per_condition/number_snps
  ][
    ,c('number_snps','number_snps_per_condition','cell_condition','fraction','pop')
  ] %>% unique()) %>% rbindlist()
x=x[
  ,fraction:=round(fraction,2)
  ]
}

immune_stimulated_rested=stimulated_vs_rested(immune_tfbs_atac)

## enrichment in these elements
open_chrom_regions=copy(immune_atacseq)[,4:=NULL] %>% unique()
open_chrom_regions=open_chrom_regions[,number_snps:=.N
                                      ][,number_snps_per_condition:=.N,by=.(cell_state)
                                      ][,cell_condition:=cell_state][,c(5:7)
                                          ][
                                            ,fraction:=number_snps_per_condition/number_snps][,pop:='Calderon et al.'] %>% unique()

immune_stimulated_rested=rbind(immune_stimulated_rested,open_chrom_regions)


test_enrichment=function(x,population){
  df=copy(x)[!pop%in%c(population,'Calderon et al.') 
             ][
               ,c(1:3,5)][
                 ,totsnps_percond:=sum(number_snps_per_condition),by=.(cell_condition)
                 ][
                   ,totnasnps:=522
                   ][
                     !pop%in%'png' 
                     ][
                       ,test:=phyper(number_snps_per_condition-1,totsnps_percond,totnasnps,number_snps,lower.tail=F)
                       
                       ]
}

nean_enrichmnet=test_enrichment(immune_stimulated_rested,'denisova')
deni_enrichmnet=test_enrichment(immune_stimulated_rested,'neandertal')

immune_stimulated_rested$pop=factor(immune_stimulated_rested$pop,levels=c('denisova','neandertal','png','Calderon et al.'))

## plot 
my_palette_cells=c('azure3', #rested
                   'lightgoldenrod3' # stimulated 
)

names(my_palette_cells)=levels(as.factor(immune_stimulated_rested$cell_condition))
colScale_cells=scale_fill_manual(name= " ",values = my_palette_cells,labels = c('Rested','Stimulated'))

cell_condition_plot=ggplot(immune_stimulated_rested,aes(x=pop,y=fraction,fill=cell_condition))+
  geom_bar(stat='identity')+colScale_cells+
  scale_x_discrete(labels=c("denisova" = "Denisova","neandertal" = "Neandertal",'png'='PNG','Calderon et al.'='Calderon et al.'))+
  ylab('\n Fraction TFBS-SNPs in open chromatin regions  \n')+
  xlab(' ')+
  theme(panel.background =element_rect(fill = 'white',size = 1, colour = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(),
        legend.title = element_text(),
        axis.text.x = element_text(),
        axis.text.y = element_text(),
        axis.title.y = element_text(hjust=0.5),
        axis.text=element_text(),
        axis.line = element_line(color = "black", size = 0, linetype = "solid"))


pdf('/home/dvespasiani/tfbs_plots/atac_immune_tfbs_snps.pdf',width = 5,height =4)
cell_condition_plot
dev.off()



