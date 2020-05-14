## CTCF location + frequencies 
library(dplyr); library(data.table)
library(magrittr); 
library(purrr)
library(GenomicRanges)
library(reshape)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/Motifbreak/')

output_dir='./CTCF_files/'

tfbs=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            drop=c('Refpvalue','Altpvalue'))[geneSymbol%like%'CTCF']
                      %>% unique() %>% dplyr::select(c('seqnames','start','end',everything()))%>% as.data.table()
      )
      x=x[c(1:3)] # dont consider naSNPs now
    
      pop_names=c('ambiguous','denisova','neandertal')
      names(x)=pop_names
      for (i in seq_along(x)){
        assign(pop_names[i],x[[i]],.GlobalEnv)}
      return(x)
}


ctcf_snps=tfbs('./Tx_and_CREs/TFBSs_disrupted_10neg5/combined/')


## get rsids
rsids=fread('../../human_genome/human_9606_b150_GRCh37p13_All_20170710_AllSNPs.gz',header = F,sep='\t',
            col.names = c('seqnames','start','rsid','ref','alt'))
rsids$seqnames=sub("^", "chr", rsids$seqnames)
rsids=rsids[!seqnames%in%c('chrMT','chrX','chrY')][,c(1:3)]


# combine
lapply(ctcf_snps,function(x)setkeyv(x,c('seqnames','start','end')))
setkeyv(rsids,c('seqnames','start'))

ctcf_snps=lapply(ctcf_snps,function(x)x=x[rsids,on=c('seqnames','start'),nomatch=0])

ctcf_snps=lapply(ctcf_snps,function(x)x=x[,delta_pwm:=abs(scoreRef-scoreAlt)] %>% group_by(`seqnames`,`start`,`end`) %>% 
              filter(`delta_pwm`==max(delta_pwm)) %>% dplyr::select(-'delta_pwm')%>% as.data.table() %>% unique())


filenames=paste0(output_dir,names(ctcf_snps),sep='')
mapply(write.table, ctcf_snps, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)



