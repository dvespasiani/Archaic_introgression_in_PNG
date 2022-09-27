## QCs output motifbreak + table SNPs
library(dplyr); library(data.table); library(magrittr)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG')

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
snps_instates=lapply(states, function(x)x=x[,c("seqnames", "start",'end','ref','alt')] %>% unique())
purrr::map2(input_snps,snps_instates,setdiff) %>% rbindlist() %>% nrow()
purrr::map2(snps_instates,input_snps,setdiff) %>% rbindlist() %>% nrow()

## read tfbs snps
read_tfbs_snps=function(x){
  file=paste0(motifbreak_path,output_snps_dir,sep='')
  file=paste0(file,x,sep='')
  
  df=as.character(list.files(file,recursive = F,full.names = T)) %>% 
    lapply(function(x)x=fread(x,sep=' ',header = T)[
      ,end:=start+1
      ][
       !dataSource %in%'HOCOMOCOv11-secondary-D'
      ]%>% setorderv(c('seqnames','start','end')) %>% 
        dplyr::select(c('seqnames','start','end','REF','ALT',everything())) %>% 
        as.data.table())
  
  pop_names=c('denisova','neandertal','png')
  names(df)=pop_names
  for (i in seq_along(df)){
    assign(pop_names[i],df[[i]],.GlobalEnv)}
  return(df)
}


jaspar_tfbs=read_tfbs_snps('jaspar2018/')
hocomoco_tfbs=read_tfbs_snps('hocomoco/')

lapply(jaspar_tfbs,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())
lapply(hocomoco_tfbs,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())

check_output_match_input=function(x,y,z){
  x=copy(x) %>% lapply(function(a)a=a[,c("seqnames", "start",'end')] %>% unique())
  y=copy(y) %>% lapply(function(a)a=a[,c("seqnames", "start",'end')] %>% unique())
  z=copy(z) %>% lapply(function(a)a=a[,c("seqnames", "start",'end')] %>% unique())
  
  df=purrr::map2(x,y,setdiff) %>%rbindlist() %>%nrow()
  df2=purrr::map2(z,y,setdiff)%>%rbindlist() %>%nrow()
  
  df3=df+df2
  if(df3%in%0){
    print('all good')
  }else{
    print('cazzo')
  }
}


check_output_match_input(jaspar_tfbs,states,hocomoco_tfbs)

### make table SNPs
## remove duplicates (i.e. SNPs that disrupt motif of same TF)
combined=purrr::map2(jaspar_tfbs,hocomoco_tfbs,rbind) %>% 
  lapply(function(x)x=x[,c('seqnames','start','end','REF','ALT','strand','geneSymbol','providerName','dataSource','alleleDiff','effect')] %>% 
           unique())

# combined=lapply(combined,function(x)x=x[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start,end,REF,ALT,geneSymbol)])

combined_dir=paste0(motifbreak_path,paste0(output_snps_dir,'combined/',sep=''),sep='')
combined_filenames=paste0(combined_dir,names(combined),sep='')
mapply(write.table,combined, file = combined_filenames,col.names = T, row.names = F, sep = " ", quote = F)

number_snps=function(x){
  df=copy(x) %>% lapply(function(y)y=y[,c('seqnames','start','end')] %>% unique()%>% nrow())
  return(df)
}

number_input=number_snps(input_snps)

number_jasp_output=number_snps(jaspar_tfbs)
number_hoco_output=number_snps(hocomoco_tfbs)
number_combined_output=number_snps(combined)


## Assign Vierstra cluster info (PS: only hocomoco v11 from secondary d are not assigned because those are absent from vierstra)
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


jaspar_cluster=tf_cluster(jaspar_tfbs)
hocomoco_cluster=tf_cluster(hocomoco_tfbs)
combined_cluster=tf_cluster(combined)

number_jasp_cluster=number_snps(jaspar_cluster)
number_hoco_cluster=number_snps(hocomoco_cluster)
number_combined_cluster=number_snps(combined_cluster)



