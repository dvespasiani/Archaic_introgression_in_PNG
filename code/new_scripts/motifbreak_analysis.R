## This script takes list of SNPs annotated within CREs and Tx states and assess their impact on known TFBSs
## I am running this script as an array (see motifbreak_commands_job_array.sh) so I/O for the script come from there
library(BSgenome)
library(motifbreakR)
library(GenomicRanges);
library(dplyr)
library(data.table)
library(magrittr)
library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocParallel)

workers = multicoreWorkers()-1

options(scipen=999)

## get command line arguments
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)


## create a for loop to iterate over the list of command line arguments
for (e in args) {
    argname <- e[1]
    if (! is.na(e[2])) {
        argval <- e[2]
        ## regular expression to delete initial \" and trailing \"
        argval <- gsub("(^\\\"|\\\"$)", "", argval)
    }
    else {
        # If arg specified without value, assume it is bool type and TRUE
        argval <- TRUE
    }
    # Infer type from last character of argname, cast val
    type <- substring(argname, nchar(argname), nchar(argname))
    if (type == "I") {
        argval <- as.integer(argval)
    }
    if (type == "N") {
        argval <- as.numeric(argval)
    }
    if (type == "L") {
        argval <- as.logical(argval)
    }
    assign(argname, argval)
    cat("Assigned", argname, "=", argval, "mode", mode(argval), "typeof", typeof(argval), "\n")
}


snps = snps.from.file(input,search.genome = BSgenome.Hsapiens.UCSC.hg19, format = "bed")

hocomoco_cores = c(
  'HOCOMOCOv11-core-A','HOCOMOCOv11-core-B','HOCOMOCOv11-core-C',
  'HOCOMOCOv11-secondary-A','HOCOMOCOv11-secondary-B','HOCOMOCOv11-secondary-C'
)


jaspar2018 = subset(MotifDb, organism=='Hsapiens' & dataSource=='jaspar2018')
hocomoco = subset(MotifDb, organism=='Hsapiens' & dataSource %in% hocomoco_cores)


motifbreak = function(file,pwmdb){
  y=motifbreakR(snpList =file, filterp = TRUE,
                pwmList = pwmdb,threshold = 1e-5,
                method = "log",bkg = c(A=0.3, C=0.2, G=0.2, T=0.3), ## use nt frequencies for non-coding regions human genome
                BPPARAM = BiocParallel::MulticoreParam(workers=workers))%>%
                as.data.table()

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

jaspar_snps_motifbreak = motifbreak(snps,jaspar2018)
hocomoco_snps_motifbreak = motifbreak(snps,hocomoco)

## write files 
write.table(jaspar_snps_motifbreak, output_jaspar ,col.names = T, row.names = F, sep = "\t", quote = F)
write.table(hocomoco_snps_motifbreak,output_hocomoco,col.names = T, row.names = F, sep = "\t", quote = F)

