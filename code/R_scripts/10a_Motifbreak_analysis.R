## motifbreak
library(BSgenome);library(motifbreakR); library(GenomicFeatures);library(GenomicRanges);
library(dplyr); library(data.table)
library(magrittr); library(IRanges); library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)

options(scipen=999)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/Motifbreak/Tx_and_CREs/')

jaspar_output_dir='./TFBSs_disrupted_10neg5/jaspar/'
hocomoco_output_dir='./TFBSs_disrupted_10neg5/hocomoco/'
encode_output_dir='./TFBSs_disrupted_10neg5/encode/'

snps=as.character(list.files('./snps_motifbreakr_format/non_splitted',recursive = F,full.names = T)) %>% 
  lapply(function(x)x=snps.from.file(file =x,
                                      search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                      format = "bed"))

encode=data('encodemotif')
hocomoco=data('hocomoco')
jaspar2018=subset(MotifDb, organism=='Hsapiens' & dataSource=='jaspar2018')


motifbreak=function(file,pwmdb){
  y=motifbreakR(snpList =file, filterp = TRUE,
              pwmList = pwmdb,threshold = 1e-5,
              method = "log",bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
              BPPARAM = BiocParallel::SerialParam())%>%as.data.table()
  }

assign_names=function(y){
  pop_names=c('ambiguous','denisova','neandertal','png')
  names(y)=pop_names
  for (i in seq_along(y)){
    assign(pop_names[i],y[[i]],.GlobalEnv)}
  return(y)
}

jaspar_snps_motifbreak=lapply(snps,function(x)x=motifbreak(x,jaspar2018))
jaspar_snps_motifbreak=assign_names(jaspar_snps_motifbreak)

encode_snps_motifbreak=lapply(snps,function(x)x=motifbreak(x,encode))
encode_snps_motifbreak=assign_names(encode_snps_motifbreak)

hocomoco_snps_motifbreak=lapply(snps,function(x)x=motifbreak(x,hocomoco))
hocomoco_snps_motifbreak=assign_names(hocomoco_snps_motifbreak)


## write files #
jaspar_filenames=paste0(jaspar_output_dir,names(jaspar_snps_motifbreak),sep='')
mapply(write.table, jaspar_snps_motifbreak, file = jaspar_filenames,col.names = T, row.names = F, sep = " ", quote = F)

hocomoco_filenames=paste0(hocomoco_output_dir,names(hocomoco_snps_motifbreak),sep='')
mapply(write.table, hocomoco_snps_motifbreak, file = hocomoco_filenames,col.names = T, row.names = F, sep = " ", quote = F)

encode_filenames=paste0(encode_output_dir,names(encode_snps_motifbreak),sep='')
mapply(write.table, encode_snps_motifbreak, file = encode_filenames,col.names = T, row.names = F, sep = " ", quote = F)




