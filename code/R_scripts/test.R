library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
motif_cluster_dir='../Annotation_and_other_files/Motifs_clusters/'
merging_keys=c('seqnames','start','end')

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)


setwd('~/Desktop/Archaic_introgression_in_PNG/')

# setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','ref','alt','chrom_state','cell_type','cell_line','all_freq'))[
        chrom_state%in%c('1_TssA','2_TssAFlnk','3_TxFlnk','4_Tx','5_TxWk',"6_EnhG","7_Enh")
        ][
          ,chrom_state:=NULL
          ] %>% unique()
    )
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set')

read_tfbs_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,end:=start+1
        ]
    )
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


combined=read_tfbs_snps(paste(motif_input_dir,'combined/',sep=''))
lapply(combined,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

## check again the input snps are the ones you wanted
combined=purrr::map2(combined,states,inner_join,by=c('seqnames','start','end','REF'='ref','ALT'='alt')) %>% 
  lapply(function(x)x=as.data.table(x))
# combined=lapply(combined,function(y)y=y[, .SD[which.min(alleleDiff)], by=.(seqnames,start)])

combined=lapply(combined,function(x)x=x[,all_freq:=round(all_freq,2)][all_freq>=0.05][,all_freq:=NULL] %>% unique())
lapply(combined,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())


test_deni=copy(combined[[1]])
# [cell_line%in%c('BCells','TCells')])
test_deni=test_deni[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='denisova'] %>% unique
test_nean=copy(combined[[2]])
# [cell_line%in%c('BCells','TCells')])
test_nean=test_nean[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='neandertal']%>% unique()

test_png=copy(combined[[3]])
# [cell_line%in%c('BCells','TCells')])
test_png=test_png[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='png']%>% unique()

test_deni=test_deni[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]
test_nean=test_nean[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]
test_png=test_png[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]

wilcox.test(test_deni$alleleDiff,test_png$alleleDiff)
wilcox.test(test_nean$alleleDiff,test_png$alleleDiff)

test=rbind(test_deni,test_nean,test_png)

ggplot(test,aes(x=pop,y=alleleDiff,fill=pop))+
  geom_violin(trim=F,scale = "width")+
  geom_boxplot(width=.4, position =  position_dodge(width = 0.4),
               outlier.size=0.2)+
  stat_compare_means(method = "wilcox",ref.group ='png')
 


test_deni_immune=copy(combined[[1]][cell_line%in%c('BCells')])
test_deni_immune=test_deni_immune[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='denisova'] %>% unique
test_nean_immune=copy(combined[[2]][cell_line%in%c('BCells')])
test_nean_immune=test_nean_immune[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='neandertal']%>% unique()

test_png_immune=copy(combined[[3]][cell_line%in%c('BCells')])
test_png_immune=test_png_immune[,c('seqnames','start','end','REF','ALT','effect','geneSymbol','providerName','alleleDiff','cell_type','cell_line')][,pop:='png']%>% unique()

test_deni_immune=test_deni_immune[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]
test_nean_immune=test_nean_immune[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]
test_png_immune=test_png_immune[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start)]


png_median_effect=median(test_png_immune$alleleDiff)
test_deni_immune_enrich=test_deni_immune[
  ,diff:=alleleDiff-png_median_effect
]
test_nean_immune_enrich=test_nean_immune[
  ,diff:=alleleDiff-png_median_effect
  ]

test_immune_enrichment=rbind(test_deni_immune_enrich,test_nean_immune_enrich)


ggplot(test_immune_enrichment,aes(x=pop,y=diff,fill=pop))+
  geom_violin(trim=F,scale = "width")+
  geom_boxplot(width=.4, position =  position_dodge(width = 0.4),
               outlier.size=0.2,notch = T)+
  stat_compare_means(method = "wilcox",method.args = list(mu=1))+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = png_median_effect)




##-------------------------------------------------------------
##    mean asnps- mean nasnps allele diff per tf
##-------------------------------------------------------------
test2=list(test_deni,test_nean,test_png)


motif=fread(paste0(motif_cluster_dir,'motifs_clusters',sep=''),sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
  ]

tf_cluster=function(x){
  df=copy(x)
  df=lapply(df,function(y)y=y[
    ,geneSymbol:=toupper(geneSymbol)
    ][
      ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub('::', '+', geneSymbol)
          ] %>%inner_join(motif,by=c('providerName'='motif_id'))%>%as.data.table())
  
  df=lapply(df,function(y)y=y[,c('seqnames','start','end','REF','ALT','alleleDiff','Cluster','Name')]%>%unique())
  return(df)
}


combined_cluster=tf_cluster(test2)


manipulate_file2=function(a){
  df=copy(a)
  df=df[,c('seqnames','start','Cluster','Name','alleleDiff')] %>% unique()
  df=df[,mean_effect_in_cluster:=mean(alleleDiff),by=.(Cluster)]
  # df=df[,c('seqnames','start'):=NULL]
  return(df)
}

combined_cluster=lapply(combined_cluster,function(x)manipulate_file2(x))

test_deni2=inner_join(combined_cluster[[1]],combined_cluster[[3]],by=c('Name','Cluster')) %>% as.data.table()
test_nean2=inner_join(combined_cluster[[2]],combined_cluster[[3]],by=c('Name','Cluster')) %>% as.data.table()

test_deni3=test_deni2[
  ,enrichment:=mean_effect_in_cluster.x/mean_effect_in_cluster.y
][
  ,c('Name','Cluster','enrichment')
] %>% unique()

test_deni3=test_deni3 %>% setorderv('enrichment',-1)

test_nean3=test_nean2[
  ,enrichment:=mean_effect_in_cluster.x/mean_effect_in_cluster.y
  ][
    ,c('Name','Cluster','enrichment')
    ] %>% unique()

test_nean3=test_nean3 %>% setorderv('enrichment',-1)


test_deni_plot=test_deni2[Cluster==161][,c('seqnames.x','start.x','alleleDiff.x','Name')][,pop:='deni'] %>% unique() %>% setnames(c('seqnames','start','alleleDiff','Name','pop'))
test_png_plot=test_deni2[Cluster==161][,c('seqnames.y','start.y','alleleDiff.y','Name')][,pop:='png'] %>% unique()%>% setnames(c('seqnames','start','alleleDiff','Name','pop'))
test_nean_plot=test_nean2[Cluster==161][,c('seqnames.y','start.y','alleleDiff.y','Name')][,pop:='nean'] %>% unique()%>% setnames(c('seqnames','start','alleleDiff','Name','pop'))

test_clust_plot=rbind(test_deni_plot,test_nean_plot,test_png_plot)


ggplot(test_clust_plot,aes(x=pop,y=alleleDiff,fill=pop))+
  geom_violin(trim=F,scale = "width")+
  geom_boxplot(width=.4, position =  position_dodge(width = 0.4),
               outlier.size=0.2)+
  stat_compare_means(method = "wilcox",ref.group ='png')



test_deni2=inner_join(combined_cluster[[1]],combined_cluster[[3]],by=c("Cluster",'Name')) %>% as.data.table()
test_deni2=test_deni2[
  ,c('seqnames.y','start.y','alleleDiff.y'):=NULL
] %>% unique()


test_deni2=test_deni2[
  ,enrichment:=alleleDiff.x/mean_effect_in_cluster.y
  ]
test_deni2=test_deni2[,c('Cluster','Name','enrichment')] %>% unique()
test_deni2=test_deni2 %>% setorderv('enrichment',-1)


ggplot(test_deni2,aes(x=enrichment))+
  geom_density()




test_deni2=copy(test_deni)
test_deni2=test_deni2

test_png2=copy(test_png)
test_png2=test_png2[
  ,mean_nasnps_alldiff:=mean(alleleDiff),by=(geneSymbol)
  ][,c('geneSymbol','mean_nasnps_alldiff','pop')] %>% unique()

test2=test_deni2[test_png2,on='geneSymbol']
test2[is.na(test2)]=0
test2=test2[
  ,c('geneSymbol','alleleDiff','mean_nasnps_alldiff','pop')
]
test2=test2[
  ,enrichment:=alleleDiff/mean_nasnps_alldiff
][
  ,c('geneSymbol','enrichment')
] %>% unique()




test2=test2[
  ,pvals:=wilcox.test(enrichment,mu=1)$p.value,by=(geneSymbol)
]
pvals=copy(test2) %>% pull(pvals)

pvalues_table=data.table(pval=pvals)
pvalues_table=pvalues_table[
  ,adj_p:=p.adjust(pval,method = 'fdr')][
    ,log10_p_adjust:=-log10(adj_p)
    ][
      ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
                                 ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
                                        ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
                                               ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
      ]

test2=cbind(test2,pvalues_table)
test2[adj_p<0.05]



  
  

ggscatter(test2, x = "mean_asnps_alldiff", y = "mean_nasnps_alldiff",
          add = "reg.line",                                 
          conf.int = TRUE,                                  
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)



ggplot(test2,aes(x=mean_asnps_alldiff))+
  geom_density()
ggplot(test2,aes(x=mean_nasnps_alldiff))+
  geom_density()
ggqqplot(test2,x='mean_nasnps_alldiff')
ggqqplot(test2,x='mean_asnps_alldiff')

## select top 5% data
n= 5
top_deni=copy(test_deni)[alleleDiff<0]
top_deni=top_deni[top_deni$alleleDiff > quantile(top_deni$alleleDiff,prob=n/100),]

top_nean=copy(test_nean)[alleleDiff<0]
top_nean=top_nean[top_nean$alleleDiff > quantile(top_nean$alleleDiff,prob=n/100),]

top_png=copy(test_png)[alleleDiff<0]
top_png=top_png[top_png$alleleDiff > quantile(top_png$alleleDiff,prob=n/100),]










motif=fread(paste0(motif_cluster_dir,'motifs_clusters',sep=''),sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
  ]

tf_cluster=function(x){
  df=copy(x)
  df=lapply(df,function(y)y=y[
    ,geneSymbol:=toupper(geneSymbol)
    ][
      ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub('::', '+', geneSymbol)
          ] %>%inner_join(motif,by=c('providerName'='motif_id'))%>%as.data.table())
  df=lapply(df,function(y)y=y[,c('seqnames','start','end','REF','ALT','effect','Cluster','Name')]%>%unique())
  return(df)
}


combined_cluster=tf_cluster(combined2)
lapply(combined_cluster,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
# lapply(combined,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())


##-------------------------------------------------------------
##    Calculate log2 fold enrichment + fisher/binomial p vals
##-------------------------------------------------------------

fisher_pvalues=function(a){
  pvals=copy(a)
  pvals=pvals[,c(1:2):=NULL]
  pvals=apply(pvals, 1, 
              function(x) {
                tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                fisher.test(tbl, alternative="g")$p.value
              })
  
  pvalues_table=data.table(pval=pvals)
  pvalues_table=pvalues_table[
    ,adj_p:=p.adjust(pval,method = 'fdr')][
      ,log10_p_adjust:=-log10(adj_p)
      ][
        ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
                                   ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
                                          ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
                                                 ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
        ]
  
  return(pvalues_table)
}

binomial_pvalues=function(a){
  binomial_test_pvalues=function(succ,fail,prob) {
    binom.test(c(succ,fail),p=prob)$p.value }
  
  a=a[
    ,pval:=mapply(binomial_test_pvalues, 
                  a$test_in_cluster,a$test_not_cluster,a$mean_ratio_cluster)
    ]
  
  pvalues=a$pval
  adj_p=p.adjust(pvalues,'fdr')
  adj_p_df=as.data.table(adj_p)
  adj_p_df=adj_p_df[
    ,log10_p_adjust:=-log10(adj_p)
    ][
      ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
                                 ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
                                        ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
                                               ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
      ]
  
  pvalues_table=cbind(a,adj_p_df)
  return(pvalues_table)
  
}

manipulate_file=function(a){
  df=copy(a)
  df=df[,c('seqnames','start','Cluster','Name')] %>% unique()
  df=df[,snps_in_cluster:=.N,by=.(Cluster)]
  tot_snps=copy(df)
  tot_snps=tot_snps[,c('seqnames','start')]%>%unique()
  tot_snps=tot_snps[,totsnps:=.N]
  
  df=inner_join(df,tot_snps,by=c('seqnames','start'))%>%setDT()
  df=df[,snps_not_cluster:=totsnps-snps_in_cluster]
  
  df=df[
    ,c('Name','Cluster','snps_in_cluster','snps_not_cluster')
    ] %>% unique()
  return(df)
}

calculate_enrichment=function(x,y){
  test=copy(x)%>%manipulate_file()
  bkgr=copy(y)%>%manipulate_file()
  df=inner_join(test,bkgr,by=c('Name','Cluster'))%>% as.data.table()
  
  df=setnames(df,old=c('snps_in_cluster.x','snps_not_cluster.x','snps_in_cluster.y','snps_not_cluster.y'),
              new=c('test_in_cluster','test_not_cluster','bkgr_in_cluster','bkgr_not_cluster'))
  df=df[
    ,bkgr_ratio:=bkgr_in_cluster/bkgr_not_cluster
    ][
      ,mean_ratio_cluster:=mean(bkgr_ratio)
      ][
        ,c('Name','Cluster','test_in_cluster','test_not_cluster',
           'bkgr_in_cluster','bkgr_not_cluster','mean_ratio_cluster')
        ] %>% unique()
  
  df_pvalues=copy(df)
  df_pvalues=fisher_pvalues(df_pvalues)
  
  df_final=cbind(df,df_pvalues)
  
  df_final=df_final[,qtl_random_ratio_cluster:=test_in_cluster/bkgr_in_cluster][
    ,mean_ratio:=mean(qtl_random_ratio_cluster)
    ][
      ,log2_fold_enrichment:=log2(qtl_random_ratio_cluster/mean_ratio)
      ][
        ,log10_numb_snps_tf:=log10(test_in_cluster) 
        ][
          ,c('Cluster','Name','log2_fold_enrichment','pval','adj_p',
             'significant_score','log10_p_adjust','log10_numb_snps_tf')
          ]%>% setorderv('log10_p_adjust',-1)
  
  df_final[mapply(is.infinite,df_final)]= NA
  df_final=df_final[
    ,maxp:=max(log10_p_adjust,na.rm=T)
    ][
      ,log10_p_adjust:=ifelse(is.na(log10_p_adjust),maxp+1,log10_p_adjust)
      ]
  
  return(df_final)  
}

deni_enrich=calculate_enrichment(combined_cluster[[1]],combined_cluster[[3]])
deni_enrich=deni_enrich[,pop:='denisova']
nean_enrich=calculate_enrichment(combined_cluster[[2]],combined_cluster[[3]])
nean_enrich=nean_enrich[,pop:='neandertal']


## table p values
export_pvalues=function(x){
  df=copy(x)[,c('Cluster','Name','log10_numb_snps_tf','pval','adj_p','log10_p_adjust','significant_score')] %>% unique()
}

denisova_pval=export_pvalues(denisova_matrix)
neandertal_pval=export_pvalues(neandertal_matrix)

pvals=list(denisova_pval,neandertal_pval)

write.xlsx(pvals,'/home/dvespasiani/Archaic_introgression_in_PNG/pvalue_tables/Supp_Table_high_freq_TFs_fisher_pvalues_corrected.xlsx')

asnps_enrichment=rbind(denisova_matrix,neandertal_matrix)

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per TF':=log10_numb_snps_tf]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
  text=ifelse(df$log10_p_adjust>=2,df$Name,'')
  
  ggplot(df,aes(x=log2_fold_enrichment,log10_p_adjust,label = text,col=log10_numb_snps_tf))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
    geom_text_repel(size = 5,color='black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines")
    )+
    gradient+
    xlab('\n Log2 fold enrichment \n')+
    ylab('-Log10 (P)')+
    xlim(-2,2)+
    facet_wrap(pop~.,ncol = 2)+
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
          axis.line = element_blank())
}


pdf('/home/dvespasiani/Archaic_introgression_in_PNG/volcano_plots/asnps_high_freq_volcano_binomial.pdf',width = 30,height=10)
tf_enrich_plot(asnps_enrichment)
dev.off()





