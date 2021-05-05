library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

options(width=150)
motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
motif_cluster_dir='../Annotation_and_other_files/Motifs_clusters/'
merging_keys=c('seqnames','start','end')

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)


plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')


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

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/')

motif=fread(paste0(motif_cluster_dir,'motifs_clusters',sep=''),sep=' ',header = T)[
  ,motif_id:=ifelse(Motif%like%'HUMAN.H11',Motif,sub(".*_", "", Motif))
  ]

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


combined=read_tfbs_snps(paste(motif_input_dir,'new_set/combined/',sep=''))
lapply(combined,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

combined=purrr::map2(combined,states,inner_join,by=c('seqnames','start','end','REF'='ref','ALT'='alt')) %>% lapply(function(x)as.data.table(x))

combined_all=copy(combined)
combined_all=lapply(combined_all,function(x)x=x[,c('cell_line','cell_type'):=NULL] %>% unique())

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
  
  df=lapply(df,function(y)y=y[,c('seqnames','start','end','REF','ALT','geneSymbol','providerName','alleleDiff','Cluster','Name')]%>%unique())
  return(df)
}


combined_all_cluster=tf_cluster(combined_all)

##-------------------------------------------------------------
##    Calculate log2 fold enrichment + fisher/binomial p vals
##-------------------------------------------------------------
## ps this is a bit nested  and convoluted. perhaps better splitting into ≠ functions but it works
fisher_pvalues=function(a){
  pvals=copy(a)
  pvals=pvals[,c(1:2):=NULL]
  pvals=apply(pvals, 1, 
              function(x) {
                tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                fisher.test(tbl)$p.value
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
  
  df_final=df_final[
    ,qtl_random_ratio_cluster:=test_in_cluster/bkgr_in_cluster
    ][
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
 
deni_enrich=calculate_enrichment(combined_all_cluster[[1]],combined_all_cluster[[3]])
deni_enrich=deni_enrich[,pop:='denisova']
nean_enrich=calculate_enrichment(combined_all_cluster[[2]],combined_all_cluster[[3]])
nean_enrich=nean_enrich[,pop:='neandertal']

asnps_enrichment=rbind(deni_enrich,nean_enrich)

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per TF':=log10_numb_snps_tf]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
  text=ifelse(!df$significant_score%in%' ',df$Name,'')
  
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


pdf(paste(plot_dir,'asnps_high_freq_fisher_volcano.pdf',sep=''),width = 10,height=6)
tf_enrich_plot(asnps_enrichment)
dev.off()

### common to high freq
combined_high=copy(combined)

combined_high=lapply(combined_high,function(x)x=x[
  ,all_freq:=round(all_freq,2)
  ][all_freq>=0.05][,c('cell_line','cell_type'):=NULL] %>% unique())

combined_high_cluster=tf_cluster(combined_high)

deni_enrich_high=calculate_enrichment(combined_high_cluster[[1]],combined_high_cluster[[3]])
deni_enrich_high=deni_enrich_high[,pop:='denisova']
nean_enrich_high=calculate_enrichment(combined_high_cluster[[2]],combined_high_cluster[[3]])
nean_enrich_high=nean_enrich_high[,pop:='neandertal']

asnps_enrichment_high=rbind(deni_enrich_high,nean_enrich_high)

## plot
pdf(paste(plot_dir,'asnps_high_freq_fisher_volcano.pdf',sep=''),width = 10,height=6)
tf_enrich_plot(asnps_enrichment_high)
dev.off()


## table p values
export_pvalues=function(x){
  df=copy(x)[,c('Cluster','Name','log2_fold_enrichment','log10_numb_snps_tf','pval','adj_p','log10_p_adjust','significant_score')] %>% unique()
}

denisova_all_pval=export_pvalues(deni_enrich)
neandertal_all_pval=export_pvalues(nean_enrich)

denisova_high_pval=export_pvalues(deni_enrich_high)
neandertal_high_pval=export_pvalues(nean_enrich_high)

pvals=list(denisova_all_pval,denisova_high_pval,neandertal_all_pval,neandertal_high_pval)
names(pvals)=c('Denisova all','Denisova common-to-high','Neanderthal all','Neanderthal common-to-high')
write.xlsx(pvals,paste(table_dir,'Supp_Table_6_TFs_fisher_pvalues.xlsx',sep=''))








