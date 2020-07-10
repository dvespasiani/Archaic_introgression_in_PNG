library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
motif_cluster_dir='./Motifbreak/Motifs_clusters/'
merging_keys=c('seqnames','start','end')

setDTthreads(8)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')


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

## remove low freq variants and see
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

## check again the input snps are the ones you wanted
combined=purrr::map2(combined,states,inner_join,by=c('seqnames','start','end','REF'='ref','ALT'='alt')) %>% 
  lapply(function(x)x=as.data.table(x)[
    ,c('cell_type','cell_line'):=NULL
  ])
combined=lapply(combined,function(x)x=x[,all_freq:=round(all_freq,2)][all_freq>=0.05][,all_freq:=NULL] %>% unique())
lapply(combined,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

## motif clusters
motif=fread(paste0(motif_cluster_dir,'motifs_clusters',sep=''),sep=' ',header = T)[
  ,geneSymbol:=ifelse(Motif%like%'H11MO',gsub("\\_HUMAN.*","",Motif),gsub("\\_MA.*","",Motif))
  ][
    ,geneSymbol:=toupper(geneSymbol)
    ][
      ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub('::', '+', geneSymbol)
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
          ] %>%inner_join(motif,by='geneSymbol')%>%as.data.table()) 
 }

combined_cluster=tf_cluster(combined)
lapply(combined_cluster,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

##-------------------------------------------------------------
##    Calculate log2 fold enrichment + fisher exact p vals
##-------------------------------------------------------------
enrichment_pvalues=function(a){
  pvals=copy(a)
  pvals=apply(pvals, 1, function(z) pval=fisher.test(matrix(as.numeric(z[3:6]),nrow=2))$p.value)
  
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



calculate_enrichment=function(x,population){
  
  manipulate_file=function(a){
  df=copy(a)
  df=lapply(df,function(b)b=b[,c('seqnames','start','end','Cluster','Name')] %>% unique())
  df=lapply(df,function(b)b=b[
    ,numbsnps_cluster:=.N,by=.(Cluster)
    ][
      ,totsnps:=.N
      ][,c('Cluster','Name','numbsnps_cluster','totsnps')])
  return(df)
  }

  
snps_in_cluster=function(a){
    df=copy(a)[,numb_asnps_in_cluster:=numbsnps_cluster.x
                     ][
                       ,numb_asnps_notin_cluster:=totsnps.x-numbsnps_cluster.x
                       ][
                         ,numb_nasnps_in_cluster:=numbsnps_cluster.y
                         ][
                           ,numb_nasnps_notin_cluster:=totsnps.y-numbsnps_cluster.y
                           ][
                             ,c('Name','Cluster','numb_asnps_in_cluster','numb_asnps_notin_cluster',
                                'numb_nasnps_in_cluster','numb_nasnps_notin_cluster')
                             ] %>% unique()
return(df)
}


final_df=function(a,b){
  
  df_final=cbind(a,b)
  
  df_final=df_final[,asnps_nasnps_ratio_cluster:=numb_asnps_in_cluster/numb_nasnps_in_cluster][
    ,mean_ratio:=mean(asnps_nasnps_ratio_cluster)
    ][
      ,log2_fold_enrichment:=log2(asnps_nasnps_ratio_cluster/mean_ratio)
      ][
        ,log10_numb_snps_tf:=log10(numb_asnps_in_cluster) 
        ][
          ,c('Cluster','Name','log2_fold_enrichment','pval','adj_p',
             'significant_score','log10_p_adjust','log10_numb_snps_tf')
          ]%>% setorderv('log10_p_adjust',-1)
  return(df_final)
}

  
  if(population=='deni'){
    
    df=copy(x)
    df=manipulate_file(df)
    df=inner_join(df[[1]],df[[3]],by=c('Name','Cluster')) %>% as.data.table()
    df=snps_in_cluster(df)
    df_pvalues=copy(df)
    df_pvalues=enrichment_pvalues(df_pvalues)
    df_final=final_df(df,df_pvalues)
    return(df_final)
  }else{
    df=copy(x)
    df=manipulate_file(df)
    df=inner_join(df[[2]],df[[3]],by=c('Name','Cluster')) %>% as.data.table()
    df=snps_in_cluster(df)
    df_pvalues=copy(df)
    df_pvalues=enrichment_pvalues(df_pvalues)
    df_final=final_df(df,df_pvalues)
    return(df_final)
  }
}

denisova_matrix=calculate_enrichment(combined_cluster,'deni')[,pop:='denisova']
neandertal_matrix=calculate_enrichment(combined_cluster,'nean')[,pop:='neandertal']

## table p values
export_pvalues=function(x){
  df=copy(x)[,c('Cluster','Name','log10_numb_snps_tf','pval','adj_p','log10_p_adjust','significant_score')] %>% unique()
}

denisova_pval=export_pvalues(denisova_matrix)
neandertal_pval=export_pvalues(neandertal_matrix)

pvals=list(denisova_pval,neandertal_pval)

write.xlsx(pvals,'/home/dvespasiani/Archaic_introgression_in_PNG/pvalue_tables/Supp_Table_high_freq_TFs_fisher_pvalues.xlsx')

asnps_enrichment=rbind(denisova_matrix,neandertal_matrix)

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per TF':=log10_numb_snps_tf]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
  text=ifelse(!df$significant_score%in% ' ',df$Name,'')
  
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


pdf('/home/dvespasiani/Archaic_introgression_in_PNG/volcano_plots/asnps_high_freq_volcano.pdf',width = 10,height=6)
tf_enrich_plot(asnps_enrichment)
dev.off()





