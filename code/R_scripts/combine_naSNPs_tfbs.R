## combine naSNPs TFBSs
library(dplyr); library(data.table); library(magrittr)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG')

motifbreak_path='./Motifbreak/Tx_and_CREs'
input_snps_dir='/snps_motifbreakr_format/non_splitted/'
output_snps_dir='/TFBSs_disrupted_10neg5/'
motif_cluster_dir='./Motifbreak/Motifs_clusters/'


read_input_snps=function(x){
  file=paste0(x,input_snps_dir,sep='')
  
  df=as.character(list.files(file,recursive = F,full.names = T)) %>% 
    lapply(function(x)x=fread(x,sep='\t',header = F,drop='V2',col.names = c('seqnames','start','name'))[
      ,start:=as.numeric(start)
      ][
        ,end:=start+1
        ][
          ,ref:=stringr::str_sub(name,start=-3,end=-3)
          ][,alt:=stringr::str_sub(name,start=-1)][,name:=NULL] %>% unique() %>% setorderv(c('seqnames','start','end')))
  
  pop_names=c('denisova','neandertal','png')
  names(df)=pop_names
  for (i in seq_along(df)){
    assign(pop_names[i],df[[i]],.GlobalEnv)}
  return(df)
}


input_snps=read_input_snps(motifbreak_path)

read_nasnps=function(x){
  
  file=paste0(motifbreak_path,output_snps_dir,sep='')
  
  df=as.character(list.files(paste0(file,x,sep=''),recursive = F,full.names = T)) %>% 
    lapply(function(a)a=fread(a,sep=' ',header = T)) %>% rbindlist() %>% setorderv(c('seqnames','start'),1)
  df=df[
    ,end:=start+1
  ]
  return(df)
}

nasnps_jaspar=read_nasnps('splitted/jaspar/')
nasnps_hocomoco=read_nasnps('splitted/hocomoco/')


check_output=function(x,y){
  x=copy(x)[,c('seqnames','start','end')] %>% unique()
  y=copy(y)[,c('seqnames','start','end')] %>% unique()
  setdiff(x,y) %>% nrow()
  
}


check_output(nasnps_jaspar,input_snps[[3]])
check_output(nasnps_hocomoco,input_snps[[3]])


dir=paste0(motifbreak_path,output_snps_dir,sep='')
write.table(nasnps_jaspar,paste0(dir,'jaspar2018/png_jaspar2018_10neg5'),sep=' ',quote = F,row.names = F,col.names = T)
write.table(nasnps_hocomoco,paste0(dir,'hocomoco/png_hocomocov11_10neg5'),sep=' ',quote = F,row.names = F,col.names = T)





