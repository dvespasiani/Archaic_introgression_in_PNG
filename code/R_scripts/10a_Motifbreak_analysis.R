## motifbreak
## this script will take days, split it by file and then run it one by one
library(BSgenome);library(motifbreakR); library(GenomicFeatures);library(GenomicRanges);
library(dplyr); library(data.table)
library(magrittr); library(IRanges); library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocParallel)

workers=multicoreWorkers()-1

options(scipen=999)
setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/Motifbreak/Tx_and_CREs/')

jaspar_output_dir='./TFBSs_disrupted_10neg5/jaspar2018/'
hocomoco_output_dir='./TFBSs_disrupted_10neg5/hocomoco/'

input_files=function(dir){
  files=list.files(dir, recursive = T,full.names=T)
  df=lapply(files,function(x)snps.from.file(x,search.genome = BSgenome.Hsapiens.UCSC.hg19, format = "bed"))
  
  file_names=as.character(list.files(dir,recursive = F,full.names = F))
  file_names=gsub("\\..*","",file_names)
  names(df)=file_names
  for (i in seq_along(df)){
    assign(file_names[i],df[[i]],.GlobalEnv)}
  return(df)
}

snps=input_files('./snps_motifbreakr_format/non_splitted/')


hocomoco_cores=c('HOCOMOCOv11-core-A','HOCOMOCOv11-core-B','HOCOMOCOv11-core-C',
                 'HOCOMOCOv11-secondary-A','HOCOMOCOv11-secondary-B','HOCOMOCOv11-secondary-C')

jaspar2018=subset(MotifDb, organism=='Hsapiens' & dataSource=='jaspar2018')
hocomoco=subset(MotifDb, organism=='Hsapiens' & dataSource%in%hocomoco_cores)


motifbreak=function(file,pwmdb){
  y=motifbreakR(snpList =file, filterp = TRUE,
                pwmList = pwmdb,threshold = 1e-5,
                method = "log",bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                BPPARAM = BiocParallel::MulticoreParam(workers=workers))%>%as.data.table()
  y$motifPos=vapply(y$motifPos, function(x) paste(x, collapse = ","), character(1L))
  y=y[
    ,motif_start:=gsub(',.*','',motifPos)
    ][
      ,motif_start:=gsub('.*-','',motif_start) %>% as.numeric()
      ][
        ,motif_end:=gsub(".*,","",motifPos)%>% as.numeric()
        ][
          ,motif_start:=start-motif_start
          ][
            ,motif_end:=end+motif_end
            ][
              ,c('end','motifPos'):=NULL
              ]
  y=y %>% dplyr::select(c(1:4,contains('motif_'),everything())) %>% as.data.table()
  return(y)
}


assign_names=function(y){
  pop_names=c('denisova','neandertal','png')
  names(y)=pop_names
  for (i in seq_along(y)){
    assign(pop_names[i],y[[i]],.GlobalEnv)}
  return(y)
}

jaspar_snps_motifbreak=lapply(snps,function(x)x=motifbreak(x,jaspar2018))
jaspar_snps_motifbreak=assign_names(jaspar_snps_motifbreak)

hocomoco_snps_motifbreak=lapply(snps,function(x)x=motifbreak(x,hocomoco))
hocomoco_snps_motifbreak=assign_names(hocomoco_snps_motifbreak)


## write files #
jaspar_filenames=paste0(jaspar_output_dir,names(jaspar_snps_motifbreak),sep='')
mapply(write.table, jaspar_snps_motifbreak, file = jaspar_filenames,col.names = T, row.names = F, sep = " ", quote = F)

hocomoco_filenames=paste0(hocomoco_output_dir,names(hocomoco_snps_motifbreak),sep='')
mapply(write.table, hocomoco_snps_motifbreak, file = hocomoco_filenames,col.names = T, row.names = F, sep = " ", quote = F)



