## Excel list of TFBS SNPs
library(data.table)
library(openxlsx)
library(magrittr)
library(dplyr)

setDTthreads(8)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/Motifbreak/')

output_dir='/home/dvespasiani/Archaic_introgression_in_PNG/tables/'


read_tfbs=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y) y=fread(y,sep=' ',header = T))
  x=lapply(x,function(y)y=y %>% dplyr::select(c('seqnames','start','end',everything())) %>% as.data.table())
  
  pop_names=c("Denisova", "Neandertal",'Papuans')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}
tfbs_snps=read_tfbs('./Tx_and_CREs/TFBSs_disrupted_10neg5/combined')
lapply(tfbs_snps,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())

## add Cluster info
motif=fread('./Motifs_clusters/motifs_clusters',sep=' ',header = T,drop='Database')[
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

tfbs_snps_cluster=tf_cluster(tfbs_snps)
lapply(tfbs_snps_cluster,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())

## add target genes and info whether they are DE or not between mentawai and koorowai
target_genes=function(targets,de){
  x=as.character(list.files(targets,recursive = F,full.names = T)) %>%
    lapply(function(y)
      y=fread(y,sep=' ',header = T,drop=c(4,5),col.names = c('seqnames','start','end','GREAT_target_gene','distanceTSS')))
  
  y=as.character(list.files(de,recursive = F,full.names = T)) %>%
    lapply(function(y)
      y=fread(y,sep=' ',header = T))
  
  z=purrr::map2(x,y,full_join,by=c('GREAT_target_gene'='gene')) %>%
    lapply(function(x)x=as.data.table(x)[
      ,DE:=ifelse(is.na(ensembl_gene_id),'no','yes')
      ][
        ,ensembl_gene_id:=NULL
        ])
  return(z)
}

tfbs_targets=target_genes('./GREAT_GO_terms/target_genes/','./GREAT_GO_terms/DE_genes_TFBS_regulated/')
lapply(tfbs_targets,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())

## get rsids
rsids=fread('../../human_genome/human_9606_b150_GRCh37p13_All_20170710_AllSNPs.gz',header = F,sep='\t',
            col.names = c('seqnames','start','rsid','ref','alt'))
rsids$seqnames=sub("^", "chr", rsids$seqnames)
rsids=rsids[!seqnames%in%c('chrMT','chrX','chrY')][,c(1:3)]


## load eQTLs
# tfbs_eqtls=as.character(list.files('./TFBS_eQTLs/',recursive = F,full.names = T)) %>%
#   lapply(function(y)
#     y=fread(y,sep=' ',header = T,select = c('seqnames','start','end','gene_seqnames','gene_start','gene_end','gene_id',
#                                             'gene_name','maf','tissue','qval',
#                                             'pval_nominal_threshold','log2_aFC','log2_aFC_lower','log2_aFC_upper'))%>% unique())
# 
# tfbs_snps=purrr::map2(tfbs_snps,tfbs_eqtls,full_join,by=c('seqnames','start','end')) %>%
#   lapply(function(x)x=as.data.table(x))


## combine with great target genes and de genes
# tfbs_snps_final=purrr::map2(tfbs_eqtls,tfbs_targets,tfbs_snps,full_join,by=c('seqnames','start','end')) %>%
#   lapply(function(x)x=x%>% dplyr::select(c(1:3,7:36,4:6)) %>%as.data.table() %>% setorderv(c('seqnames','start','end'),1))
# 
# tfbs_snps_final=lapply(tfbs_snps_final,function(x)x=inner_join(x,rsids,by=c('seqnames','start')) %>% as.data.table())
tfbs_snps_final=purrr::map2(tfbs_snps_cluster,tfbs_targets,full_join,by=c('seqnames','start','end'))%>%
  lapply(function(x)x=x %>% as.data.table()%>% setorderv(c('seqnames','start','end'),1))

tfbs_with_rsids=copy(tfbs_snps_final) %>% lapply(function(x)x=x[rsids,on=c('seqnames','start'),nomatch=0])

tfbs_without_rsids=copy(tfbs_snps_final)
tfbs_without_rsids=purrr::map2(tfbs_without_rsids,tfbs_with_rsids,anti_join,by=c('seqnames','start','end')) %>% 
  lapply(function(x)x=as.data.table(x)[,'rsid':=NA]) 

tfbs_snps_final=purrr::map2(tfbs_with_rsids,tfbs_without_rsids,rbind)

assign_names=function(y){
  pop_names=c('Supp_file_Denisova_TFBS_aSNPs','Supp_file_Neanderta_TFBS_aSNPs','Supp_file_PNG_TFBS_naSNPs')
  names(y)=pop_names
  for (i in seq_along(y)){
    assign(pop_names[i],y[[i]],.GlobalEnv)}
  return(y)
}


tfbs_snps_final=assign_names(tfbs_snps_final)

tfbs_filenames=paste0(output_dir,names(tfbs_snps_final),sep='')
mapply(write.table, tfbs_snps_final, file = tfbs_filenames,col.names = T, row.names = F, sep = "\t", quote = F)


