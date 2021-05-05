## script used to label aSNPs and naSNPs 
## files used here dont have Baining samples

library(dplyr);
library(magrittr);
library(data.table)
library(ggplot2)
library(openxlsx)

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
out_dir = './filtered_files/'
plot_dir = './Results/Plots/QCs/snps_qcs/'
table_output = './Results/Tables/'

input_dir = './Original_files'
options(width = 150)

range_keys = c('seqnames','start','end')


##------------
## read files
##------------
## for each snp look whether it is present within archaic haps 

all_snps = list.files(input_dir,recursive =F,full.names = T,pattern='SNPstats') %>% 
    lapply(function(y)
      fread(y,sep = '\t',header = T)[
        POP_ARCH_REF + POP_ARCH_ALT>0,ancestry := 'archaic' 
        ][
          is.na(ancestry), ancestry := 'non_archaic'
          ]%>%setnames(old=c('CHR','FROM','TO'),new = range_keys)
)
names(all_snps) = c('denisova_file','neanderthal_file')

## separate naSNPs from aSNPs 
## and get the naSNPs that are in common between the two files

naSNPs = copy(all_snps)%>% lapply(function(x)x=x[ancestry=='non_archaic'][,c(11:14):= NULL])
naSNPs = naSNPs[[1]][naSNPs[[2]],on=c(range_keys,'REF',"ALT",'ANC','POP_ARCH_REF','POP_ARCH_ALT','POP_NOTARCH_REF','POP_NOTARCH_ALT','ancestry'),nomatch=0]

d_aSNPs = copy(all_snps[[1]])[!naSNPs,on=range_keys][ancestry != 'non_archaic'] ## this removes those instances called as aSNPs in the other file
n_aSNPs = copy(all_snps[[2]])[!naSNPs,on=range_keys][ancestry != 'non_archaic']


## create single list with all sorted SNPs
## remove SNPs with unknown ANC 
## calculate MIAF/DAF and remove rare variants
## determine whether the main alleles are derived or ancestral
sorted_SNPs = list(d_aSNPs,n_aSNPs,naSNPs)%>%lapply(function(x)x=x[
  !ANC=='-1'
  ][
    ,MAF := ifelse(ancestry=='archaic', ## MIAF for aSNPs
      ifelse(POP_ARCH_REF > POP_ARCH_ALT,
        POP_ARCH_REF / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
        POP_ARCH_ALT / (POP_ARCH_REF + POP_ARCH_ALT + POP_NOTARCH_REF + POP_NOTARCH_ALT)
      ), ## DAF for naSNPs
      ifelse(ANC == 0, 
      POP_NOTARCH_ALT / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
      POP_NOTARCH_REF / (`POP_ARCH_REF`+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT))
      )
  ][
    ,state_allele := ifelse(ancestry=='archaic', ## label the state of the main introgressed archaic allele
      ifelse(
        (ANC==0 & POP_ARCH_REF>POP_ARCH_ALT)|(ANC==1 & POP_ARCH_ALT>POP_ARCH_REF),'ancestral',
        ifelse(POP_ARCH_REF==POP_ARCH_ALT,'not_determined','derived')
      ),
      'derived' ## all naSNPs that will be considered are derived. SNPs fixed for the ancestral allele are removed
    )
  ][
    MAF >= 0.05 ## keep only variants having MIAF/DAF >= 0.05
  ]
)
names(sorted_SNPs) = c('denisova','neanderthal','humans')


## remove SNPs segregating within africans
## read files with the allele frequencies for the 1KG African SNPs
## then remove those variants segregating in africans at freq > 0.005
afr_1kg_snp_freq  = fread('../Annotation_and_other_files/human_genome/1kg_snp_modified.gz',sep=' ',header=T,select=c(1:5))
colnames(afr_1kg_snp_freq)[1:2] = range_keys[-3]

snps_notin_afr = lapply(sorted_SNPs,function(x) x=x[afr_1kg_snp_freq,on=c(range_keys[-3],'REF','ALT'),nomatch=0][AFR_AF<=0.005][,AFR_AF:=NULL])

snps_notin_1kg = lapply(sorted_SNPs,function(x)x=x[!afr_1kg_snp_freq,on=c(range_keys[-3],'REF','ALT')])

## Now from the filtered set of SNPs
## remove aSNPs that are in common btw arch and non-arch haplotypes
afr_filtered_snps = Map(rbind,snps_notin_afr,snps_notin_1kg)%>%
lapply(function(x)
  x=x[
    ,snp_distribution_across_haps := ifelse(
      ancestry=='non_archaic',1,
      ifelse(POP_ARCH_REF > POP_ARCH_ALT,
      (POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT)) - (POP_NOTARCH_REF/(POP_NOTARCH_REF+POP_NOTARCH_ALT)),
      (POP_ARCH_ALT/(POP_ARCH_REF+POP_ARCH_ALT)) - (POP_NOTARCH_ALT/(POP_NOTARCH_REF+POP_NOTARCH_ALT))
      )
    )
  ][
    snp_distribution_across_haps >= 0.25
  ][
    ,snp_distribution_across_haps:=NULL
  ]
)

## remove ambiguous aSNPs, i.e. those in common between deni and nean
ambiguous_asnps = afr_filtered_snps[[1]][
  afr_filtered_snps[[2]],on=c(range_keys,'REF',"ALT",'ANC','ancestry','DENI_REF','DENI_ALT','NEAN_REF','NEAN_ALT'),nomatch=0
  ][
    , MAF := ifelse(MAF > i.MAF, MAF,i.MAF)
    ][
      ,POP_ARCH_REF:=ifelse(MAF >i.MAF, POP_ARCH_REF,i.POP_ARCH_REF)
      ][
        ,POP_ARCH_ALT:=ifelse(MAF >i.MAF, POP_ARCH_ALT,i.POP_ARCH_ALT)
        ][
          ,POP_NOTARCH_REF:=ifelse(MAF >i.MAF, POP_NOTARCH_REF,i.POP_NOTARCH_REF)
          ][
            ,POP_NOTARCH_ALT:=ifelse(MAF >i.MAF, POP_NOTARCH_ALT,i.POP_NOTARCH_ALT)
            ][
              ,main_introgr_snp := ifelse(POP_ARCH_REF>POP_ARCH_ALT,'ref','alt')
              ][
                ,deni_state :=ifelse(DENI_REF>DENI_ALT,'ref','alt')
                ][
                  ,nean_state :=ifelse(NEAN_REF>NEAN_ALT,'ref','alt')
                  ][
                    ,match_genome := ifelse(
                      (main_introgr_snp==deni_state & main_introgr_snp==nean_state) | (main_introgr_snp!=deni_state & main_introgr_snp!=nean_state),
                      'ambiguous',ifelse((main_introgr_snp==deni_state & main_introgr_snp!=nean_state),'denisova','neanderthal'))
                      ]%>% dplyr::select(-c(contains('i.'))
)


aSNPs = copy(afr_filtered_snps[c(1:2)])%>%lapply(function(x)x=x[!ambiguous_asnps,on=range_keys])

get_archaic_snps = function(x,y,genome){
  archaic_snps = rbind(copy(x),copy(y)[match_genome == genome][,c(18:21):=NULL])[
    ,introgressed_allele:=ifelse(POP_ARCH_REF>POP_ARCH_ALT,'ref','alt')
    ][
      ,c(..range_keys,"REF",'ALT','MAF','introgressed_allele','state_allele')
    ]
  return(archaic_snps)
}

denisova_aSNPs = get_archaic_snps(aSNPs[[1]],ambiguous_asnps,'denisova')
neanderthal_aSNPs = get_archaic_snps(aSNPs[[2]],ambiguous_asnps,'neanderthal')

nasnps = copy(afr_filtered_snps[[3]])[,derived_allele:=ifelse(ANC==0,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','derived_allele','state_allele')]

final_set = list(denisova_aSNPs,neanderthal_aSNPs,nasnps)
names(final_set) = c('denisova','neanderthal','modern_humans')


filenames = paste0(out_dir,names(final_set),sep='.txt')
mapply(write.table,final_set, file = filenames,col.names = T, row.names = F, sep = "\t", quote = F)

