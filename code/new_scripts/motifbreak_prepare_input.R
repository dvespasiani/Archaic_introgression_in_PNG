#### Motifbreak for all SNPs
library(BSgenome)
library(motifbreakR)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(magrittr)
library(IRanges)
library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)

options(scipen=999)

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

tx_cres_states = c('1_TssA','2_TssAFlnk','3_TxFlnk',"4_Tx","5_TxWk",'6_EnhG','7_Enh')
input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
output_dir = './Motifbreak/input_files/'

cres_snps = list.files(input_dir,recursive = F,full.names = T,pattern='new_*') %>%
    lapply(function(y)fread(y,sep='\t',header=T)
    [
      chrom_state %in% tx_cres_states
      ][
        ,c('seqnames','start','end','REF','ALT')
        ] %>% unique()
)
names(cres_snps) = c('denisova','modern_humans','neanderthal')

motifbreak_format = copy(cres_snps)%>%
  lapply(function(x)
  x=x[
    ,Start:= start-1
    ][
      ,End:=start
      ][
        ,name:=paste(seqnames,End,REF,ALT,sep=':')
        ][
          ,c('seqnames','Start','End','name')] %>% unique()
)

## write files 
filenames = paste0(output_dir,names(motifbreak_format),'.bed',sep='')
mapply(write.table,motifbreak_format, file = filenames,col.names = F, row.names = F, sep = "\t", quote = F)

