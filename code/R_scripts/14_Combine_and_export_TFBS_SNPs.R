## Excel list of TFBS SNPs
library(data.table)
library(openxlsx)
library(magrittr)
library(dplyr)

options(width = 150)
numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

output_dir='./Results/Tables/'


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
tfbs_snps=read_tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/new_set/combined')
lapply(tfbs_snps,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())

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

tfbs_targets=target_genes('./Motifbreak/GREAT_GO_terms/target_genes/','./Motifbreak/GREAT_GO_terms/DE_genes_TFBS_regulated/')

## get rsids
rsids=fread('../Annotation_and_other_files/human_genome/human_9606_b150_GRCh37p13_All_20170710_AllSNPs.gz',header = F,sep='\t',
            col.names = c('seqnames','start','rsid','ref','alt'))
rsids$seqnames=sub("^", "chr", rsids$seqnames)
rsids=rsids[!seqnames%in%c('chrMT','chrX','chrY')][,c(1:3)]

## combine files 
tfbs_snps_final=purrr::map2(tfbs_snps,tfbs_targets,full_join,by=c('seqnames','start','end'))%>%
  lapply(function(x)x=x %>% as.data.table()%>% setorderv(c('seqnames','start','end'),1))

tfbs_with_rsids=copy(tfbs_snps_final) %>% lapply(function(x)x=x[rsids,on=c('seqnames','start'),nomatch=0])

tfbs_without_rsids=copy(tfbs_snps_final)
tfbs_without_rsids=purrr::map2(tfbs_without_rsids,tfbs_with_rsids,anti_join,by=c('seqnames','start','end')) %>% 
  lapply(function(x)x=as.data.table(x)[,'rsid':=NA]) 

tfbs_snps_final=purrr::map2(tfbs_with_rsids,tfbs_without_rsids,rbind)

names(tfbs_snps_final)=c('Supp_file_Denisovan_TFBS_aSNPs','Supp_file_Neandertal_TFBS_aSNPs','Supp_file_PNG_TFBS_naSNPs')

tfbs_filenames=paste0(output_dir,names(tfbs_snps_final),sep='')
mapply(write.table, tfbs_snps_final, file = tfbs_filenames,col.names = T, row.names = F, sep = "\t", quote = F)


