#### Motifbreak for all SNPs
library(BSgenome);library(motifbreakR); library(GenomicFeatures);library(GenomicRanges);
library(dplyr); library(data.table)
library(magrittr); library(IRanges); library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)

options(scipen=999)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')
snps_output='./Motifbreak/Tx_and_CREs/snps_motifbreakr_format/non_splitted/'

tx_cres_states=c('1_TssA','2_TssAFlnk','3_TxFlnk',"4_Tx","5_TxWk",'6_EnhG','7_Enh')

read_cres=function(x){
  df=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
           fread(y,sep=' ',header=T)[
             chrom_state%in%tx_cres_states
           ][,c('seqnames','start','end','ref','alt')] %>% unique()
    )
  
  pop_names=as.character(list.files(x,full.names = F,recursive = F))
  names(df)=c('denisova','neandertal','png')
  for (i in seq_along(df)){
    assign(pop_names[i],df[[i]],.GlobalEnv)}
  return(df)
  }

snps_cres=read_cres('./Chromatin_states/SNPs_chromHMM_annotated/new_set/')

# create motifbreakR format which is chr:start:REF:ALT
# remember to substract one bp from start because otherwise it understands un cazzo 
motifbreak_format= function(x,name_column){
  x=as.data.table(x)
  x=x[,Start:= start-1
      ][
        ,End:=start
      ][
        ,name:=paste(seqnames,End,ref,alt,sep=':')
      ]
  name_column=copy(x)
  name_column=name_column[
        ,c('name')
      ]
  x=do.call(cbind,list(x,name_column))
  x=x[,c('seqnames','Start','End','name')] %>% unique()

}

snps_motifbreak=lapply(snps_cres,function(x)motifbreak_format(x))

#write.files '\t'
filenames=paste0(snps_output,names(snps_motifbreak),sep='')
mapply(write.table,snps_motifbreak, file = filenames,col.names = F, row.names = F, sep = "\t", quote = F)
