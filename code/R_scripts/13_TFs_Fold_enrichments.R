library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

setDTthreads(10)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

merging_keys=c('seqnames','start','end')
motif_cluster_dir='./Motifbreak/Motifs_clusters/'

gtex_expr=fread('./GTEx/GTEx_Analysis_v8_gene_median_tmp/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',sep='\t',header = T)
gene_ids=gtex_expr[,c(1:2)] 

## get motif clusters

## Add vierstra clusters
clusters=function(clusters,cluster_names){
  df=fread(paste0(motif_cluster_dir,clusters,sep=''),sep=' ',header = T,select = c('Cluster','Motif'))
  
  df_names=fread(paste0(motif_cluster_dir,'motif_cluster_names',sep=''),sep='\t',header = T,select = c('Name','Cluster','Seed_motif'))[
    ,motif_name:=gsub('_.*','',Seed_motif)
    ][
      ,Seed_motif:=NULL]
  
  df_combined=inner_join(df,df_names,by='Cluster') %>% as.data.table()
  df_combined=df_combined[,symbol:=toupper(gsub('_.*','',Motif))]
  
  df_combined1=copy(df_combined)[,c(1:4)]
  df_combined2=copy(df_combined)[,4:=NULL]
  colnames(df_combined2)[4]='motif_name'
  
  df_combined_final=rbind(df_combined1,df_combined2) %>% unique()
  df_combined_final=df_combined_final[!motif_name%like%'mouse' ][!motif_name%like%'MOUSE']
  
  df_combined_final=semi_join(df_combined_final,gene_ids,by=c('motif_name'='Description')) %>%as.data.table() %>% unique()
  return(df_combined_final)
}


motif_clusters=clusters('motifs_vierstra','motif_cluster_names')
# motifs_clusters=fread('./Motifbreak/Motifs_clusters/motifs_vierstra',sep=' ',header = T,select = c('Cluster','geneSymbol'))
# motifs_clusters=semi_join(motifs_clusters,gene_ids,by=c('geneSymbol'='Description')) %>%as.data.table()
# 
motifs_cluster_names=fread('./Motifbreak/Motifs_clusters/motif_cluster_names',sep='\t',header = T,select = c('Name','Cluster','Seed_motif'))
motifs_cluster_names=motifs_cluster_names[!Seed_motif%like%'mouse'][!Seed_motif%like%'MOUSE'][,Seed_motif:=gsub('_.*','',Seed_motif)][!Seed_motif%like%'mouse']


tfbs=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            select =c('seqnames','start','end','REF','ALT','geneSymbol','providerName','effect')) %>% unique())
  x=lapply(x,function(y)
    y=y[
      ,duplicated_effect:=.N,by=.(seqnames,start,end,geneSymbol)
      ]%>% group_by(seqnames,start,end,geneSymbol) %>%
      filter(duplicated_effect==max(duplicated_effect)) %>% as.data.table()
  )
  x=x[c(2:4)] # remove ambiguous like this for now
  x=lapply(x,function(y)y=semi_join(y,gene_ids,by=c('geneSymbol'='Description')) %>% ## remove TFs absent from gtex (only D/OBOXs)
                         as.data.table()) 
  
  x=lapply(x,function(y)y=y[!providerName%like%'disc'][,c('providerName','effect','duplicated_effect'):=NULL][
    ,geneSymbol:=toupper(geneSymbol)
    ][
      , geneSymbol:=gsub('::', '+', geneSymbol)
      ] %>% unique())
  
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

snps_tfbs=tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/combined/')


snps_tfbs=lapply(snps_tfbs,function(x)inner_join(x,motif_clusters,by=c('geneSymbol'='motif_name'))%>% as.data.table() %>% unique())

## count number SNPs per pop per tf
snp_matrix=function(x,pop){
  snp_count=function(x){
    df=copy(x)
   
    df=lapply(df,function(x)x=x[,c('seqnames','start','end','Cluster')] %>% unique())
    df=lapply(df,function(y)y=unique(y)[,numbsnps:=.N,by=.(Cluster)][,totsnps:=.N][,c('Cluster','numbsnps','totsnps')] %>% unique())
    df=Map(mutate,df,'pop'=names(df))
    
    df_final=full_join(df[[1]],df[[2]],by='Cluster') %>%as.data.table() %>% na.omit()
    
    df_final=df_final[,tot_asnps_tf:=numbsnps.x
                      ][
                        ,tot_nasnps_tf:=numbsnps.y
                        ][
                          ,tot_asnps:=totsnps.x
                          ][
                            ,tot_nasnps:=totsnps.y
                            ][
                              ,tot_snps_tf:=numbsnps.x+numbsnps.y
                              ][
                                ,totsnps:=sum(numbsnps.x)+sum(numbsnps.y)
                                ][
                                  ,log2_fold_enrichment:=round(log2(((tot_asnps_tf/tot_asnps)/(tot_nasnps_tf/tot_nasnps))),2)
                                  ][
                                    ,c('Cluster','log2_fold_enrichment','tot_asnps_tf','tot_asnps','tot_nasnps','tot_nasnps_tf','tot_snps_tf')
                                    ]%>% unique()
    
    df_final=df_final[
      ,pval:=phyper(tot_asnps_tf-1,tot_snps_tf,tot_nasnps,tot_asnps,lower.tail = F) # this calculates the p value of having a >= numb aSNPs targeting a tf given a total number of aSNPs, naSNPs
      ][
        ,adj_p:=p.adjust(pval, method="bonferroni")
        ][
          ,log10_p_adjust:=-log10(adj_p)
          ][
            ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
                                       ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
                                              ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
                                                     ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
            ][
              ,log10_tot_asnps_tf:=log10(tot_asnps_tf)
              ][
                ,c('Cluster','log2_fold_enrichment','pval','adj_p','significant_score','log10_p_adjust','log10_tot_asnps_tf')
                ]%>% setorderv('log10_p_adjust',-1)
    
    
    df_final=inner_join(df_final,motifs_cluster_names,by='Cluster') %>% as.data.table()
    
    
  }
  
  if(pop=='deni'){
    y=copy(x)
    y=y[c(1,3)]
    y=snp_count(y)
    return(y)
  }else{
    z=copy(x)
    z=z[c(2,3)]
    z=snp_count(z)
    return(z)
  }
  
}


denisova_matrix=snp_matrix(snps_tfbs,'deni')
neandertal_matrix=snp_matrix(snps_tfbs,'nean')


## table p values
export_pvalues=function(x){
df=copy(x)[,c(1,8,9,2:5)] %>% unique()
}

denisova_pval=export_pvalues(denisova_matrix)
neandertal_pval=export_pvalues(neandertal_matrix)

pvals=list(denisova_pval,neandertal_pval)

write.xlsx(pvals,'/home/dvespasiani/pvalue_tables/Supp_Table_TFs_hypergeom_pvalues.xlsx')

## plot
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number aSNPs per TF':=log10_tot_asnps_tf]
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
  text=ifelse(df$log2_fold_enrichment>0 & df$log10_p_adjust>10,df$Name,'')
  
  ggplot(df,aes(x=log2_fold_enrichment,log10_p_adjust,label = text,col=log10_tot_asnps_tf))+
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
    theme(
      legend.position = "bottom",
      legend.key = element_rect(fill = "white", colour = "black"),
      panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
     axis.line = element_line(color = "black", size = 0, linetype = "solid"))
}


pdf('/home/dvespasiani/volcano_plots/deni_volcano.pdf',width = 7,height=7)
tf_enrich_plot(denisova_matrix)
dev.off()

pdf('/home/dvespasiani/volcano_plots/nean_volcano.pdf',width = 7,height=7)
tf_enrich_plot(neandertal_matrix)
dev.off()




