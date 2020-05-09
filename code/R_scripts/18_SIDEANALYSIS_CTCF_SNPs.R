## CTCF location + frequencies 
library(dplyr); library(data.table)
library(magrittr); 
library(purrr)
library(GenomicRanges)
library(reshape)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

tfbs=function(x){
  as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            drop=c('Refpvalue','Altpvalue','dataSource'))[
              ,start:=snpPos
              ][
                ,end:=start+1
                ][
                  ,c('snpPos'):=NULL
                  ][geneSymbol%in%c('CTCFL','CTCF')]
                      %>% unique() 
    )
}

assign_pop_names=function(x){
  pop_names=c('ambiguous_ctcf','denisova_ctcf','neandertal_ctcf','png_ctcf')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

ctcf_snps=tfbs('./Motifbreak/Temporary_tx_cres_combined/')
ctcf_snps=assign_pop_names(ctcf_snps)

# read original files
read_snps=function(x){
  as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,select = c('seqnames','start','end','ANC','all_frequency','freq_range','CADD_PHRED','CADD_RAW')
      )
    )
}

ooa_snps=read_snps('./OoA_snps')

# combine
keys=c('seqnames','start','end')
lapply(ctcf_snps,function(x)setkeyv(x,keys))
lapply(ooa_snps,function(x)setkeyv(x,keys))

ctcf_snps=purrr::map2(ctcf_snps,ooa_snps,merge,by=keys)
ctcf_snps=lapply(ctcf_snps,function(x)x=x[,c(1:3,5:7,20,21:24,9:19)] %>% unique())

write_list=function(x){
lapply(names(x),function(y, x)
  write.table(x[[y]],paste(y, "", sep = ""),row.names=FALSE, sep="\t", quote=F),x)
}


setwd('/home/dvespasiani/')
write_list(ctcf_snps)

