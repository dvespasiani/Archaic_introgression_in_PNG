#### Motifbreak for all SNPs
library(BSgenome);library(motifbreakR); library(GenomicFeatures);library(GenomicRanges);
library(dplyr); library(data.table)
library(magrittr); library(IRanges); library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)

options(scipen=999)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/Chromatin_states/SNPs_chromatin_states_without_densities')

tx_cres_states=c('1_TssA','2_TssAFlnk','3_TxFlnk',"4_Tx","5_TxWk",'6_EnhG','7_Enh')

read_cres=function(x){
  as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
           fread(y,sep=' ',header=T)[
             chrom_state%in%tx_cres_states
           ]
    )
  }

assign_names=function(x,list_df){
  pop_names=as.character(list.files(x,full.names = F,recursive = F))
  # pop_names=gsub("_chromatin_states.*","\\1",pop_names)
  names(list_df)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(list_df)
}

asnps_cres_states=read_cres('./aSNPs')
asnps_cres_states=assign_names('./aSNPs',asnps_cres_states)
asnps_cres_states=Map(mutate, asnps_cres_states,"pop"=names(asnps_cres_states))

nasnps_cres_states=fread('./naSNPs/non_archaics_png',sep=' ',header=T)[
  chrom_state%in%tx_cres_states
  ]

nasnps_cres_states=nasnps_cres_states %>% split(as.factor(nasnps_cres_states$freq_range))

assign_nasnps_names=function(x){
  nasnps_names=c('nasnps_high','nasnps_low')
      names(x)=nasnps_names
  for (i in seq_along(x)){
    assign(nasnps_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

nasnps_cres_states=assign_nasnps_names(nasnps_cres_states)
nasnps_cres_states=Map(mutate, nasnps_cres_states,"pop"=names(nasnps_cres_states))

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

asnps_motifbreak=lapply(asnps_cres_states,function(x)motifbreak_format(x))
asnps_motifbreak=assign_names('./aSNPs',asnps_motifbreak)

nasnps_motifbreak=lapply(nasnps_cres_states,function(x) motifbreak_format(x))
nasnps_motifbreak=assign_nasnps_names(nasnps_motifbreak)

#write.files '\t'
setwd('../../Motifbreak/snps_motifbreakr_format')

lapply(names(asnps_motifbreak),function(x, asnps_motifbreak)
  write.table(asnps_motifbreak[[x]], paste(x, "", sep = ""),row.names = F,col.names = F, sep="\t", quote=F),asnps_motifbreak)

lapply(names(nasnps_motifbreak),function(x, nasnps_motifbreak)
  write.table(nasnps_motifbreak[[x]], paste(x, "", sep = ""),row.names = F,col.names = F, sep="\t", quote=F),nasnps_motifbreak)
