library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
motif_cluster_dir='./Motifbreak/Motifs_clusters/'


setDTthreads(8)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')


gtex_expr=fread('./GTEx/GTEx_Analysis_v8_gene_median_tmp/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',sep='\t',header = T)
gene_ids=gtex_expr[,c(1:2)]


read_tfbs_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,c('start', 'end') :=NULL][
          ,start := snpPos][
            ,end := start +1][
              ,snpPos := NULL]
    )
  x=x[c(2:4)]# remove ambig
  x=lapply(x,function(y)y=semi_join(y,gene_ids,by=c('geneSymbol'='Description')) %>% ## remove TFs absent from gtex (only D/OBOXs)
             as.data.table()) 
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


hocomoco=read_tfbs_snps(paste(motif_input_dir,'hocomoco',sep=''))
jaspar=read_tfbs_snps(paste(motif_input_dir,'jaspar',sep=''))
encode=read_tfbs_snps(paste(motif_input_dir,'encode',sep='')) %>% lapply(function(x)x=copy(x)[providerName%like%'_known_'])


## motif clusters

motif=fread(paste0(motif_cluster_dir,'motifs',sep=''),sep='\t',header = T)[
  ,motif_names:=toupper(gsub('_.*','',Motif))
  ][!motif_names%like%'MOUSE'][!motif_names%like%'mouse']

motif_id=fread(paste0(motif_cluster_dir,'motifs_ids',sep=''),sep='\t',header = T)

motifs=inner_join(motif,motif_id,by='Cluster') %>% as.data.table()


tf_cluster=function(x,db){
  # df=rbind(jaspar_snps_motifbreak,encode_snps_motifbreak,hocomoco_snps_motifbreak)
  # df$motifPos=vapply(df$motifPos, function(x) paste(x, collapse = ","), character(1L))
  df=copy(x)
  df=lapply(df,function(y)y=y[
      ,geneSymbol:=toupper(geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
          ][
            ,geneSymbol:=gsub('::', '+', geneSymbol)
            ][
              ,providerName:=gsub('::', '+', providerName)
              ][!providerName%like%'disc'][
                ,providerName:=gsub('_.*','',providerName)
                ] %>% unique()  )
  
  combine_df=function(x){
    df1=copy(x)[,geneSymbol:=NULL]
    df2=copy(x)[,providerName:=NULL]
    names(df2)[names(df2) == 'geneSymbol'] ='providerName'
    df_combined=rbind(df1,df2)
    return(df_combined)
  }
  
  # df$geneSymbol= gsub('::', '+',  df$geneSymbol)
  # df=df[
  #   ,motif_start:=gsub(',.*','',motifPos)
  #   ][
  #     ,motif_start:=gsub('.*-','',motif_start) %>% as.numeric()
  #     ][
  #       ,motif_end:=gsub(".*,","",motifPos)%>% as.numeric()
  #       ][
  #         ,motif_start:=start-motif_start
  #         ][
  #           ,motif_end:=end+motif_end
  #           ][
  #             ,c('end','motifPos'):=NULL
  #             ][!providerName%like%'disc'] %>% unique()
  if(db=='jaspar'){
    df_jaspar=lapply(df,function(y)y=inner_join(y,motifs,by=c('geneSymbol'='motif_names'))%>%as.data.table())
    df_jaspar=lapply(df_jaspar,function(y)y=y[,providerName:=NULL][,providerName:=geneSymbol][,geneSymbol:=NULL] %>% 
                       dplyr::select(c(1:9,'dataSource','providerName',everything())) %>%
                       as.data.table()%>% setorderv(c('seqnames','start'),1))
    return(df_jaspar)
  }else{
    
    df=lapply(df,function(x)x=combine_df(x))
    
    df_final=lapply(df,function(y)y=inner_join(y,motifs,by=c('providerName'='motif_names')) %>% dplyr::select(c(1,2,3,everything())) %>%
                      as.data.table()%>% setorderv(c('seqnames','start'),1))
    return(df_final)
  }
}

hocomoco_cluster=tf_cluster(hocomoco,'hocomoco')
lapply(hocomoco_cluster,function(x)x[,c(1:3)] %>% unique() %>% nrow())
lapply(hocomoco,function(x)x[,c(1:3)] %>% unique() %>% nrow())

jaspar_cluster=tf_cluster(jaspar,'jaspar')
lapply(jaspar_cluster,function(x)x[,c(1:3)] %>% unique() %>% nrow())
lapply(jaspar,function(x)x[,c(1:3)] %>% unique() %>% nrow())

encode_cluster=tf_cluster(encode,'encode')
lapply(encode_cluster,function(x)x[,c(1:3)] %>% unique() %>% nrow())
lapply(encode,function(x)x[,c(1:3)] %>% unique() %>% nrow())


## combine files
combine=function(x,y,z){
  df=purrr::map2(x,y,rbind) %>% lapply(function(x)setDT(x))
  df=purrr::map2(df,z,rbind) %>% lapply(function(x)setDT(x))
  df=lapply(df,function(x)x=x[,'Motif':=NULL] %>%
              unique() %>%dplyr::select(c('seqnames','start','end','REF','ALT',everything())) %>% as.data.table())
}

combined_tfbs=combine(encode_cluster,jaspar_cluster,hocomoco_cluster)


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
    
    
    df_final=inner_join(df_final,motif_id,by='Cluster') %>% as.data.table()
    
    
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


denisova_matrix=snp_matrix(combined_tfbs,'deni')
neandertal_matrix=snp_matrix(combined_tfbs,'nean')


## table p values
export_pvalues=function(x){
  df=copy(x)[,c('Cluster','Name','log10_tot_asnps_tf','log2_fold_enrichment','pval','adj_p','log10_p_adjust','significant_score')] %>% unique()
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




