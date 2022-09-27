## script used to label aSNPs and naSNPs 
## files used here dont have Baining samples

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(openxlsx)

setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
out_dir = './filtered_files/'
plot_dir = './Results/Plots/QCs/'
table_output = './Results/Tables/'
input_dir = './Original_files'
options(width = 150)

source('./scripts/reusable_functions.R')
##------------
## read files
##------------
## read b-statistic values 
bstatistic <- fread('../Annotation_and_other_files/bkgd_values/bkgd_hg19.bed.gz',header=F,drop='V4',col.names=c(range_keys,'bvalue'))
setkeyv(bstatistic,range_keys)

## for each snp look whether it is present within archaic haps 
all_snps = list.files(input_dir,recursive =F,full.names = T,pattern='UVBaining') %>% 
    lapply(function(y){
      y <-fread(y,sep = '\t',header = T)[
        POP_ARCH_REF + POP_ARCH_ALT>0,ancestry := 'archaic' 
        ][
          is.na(ancestry), ancestry := 'non_archaic'
          ]%>%setnames(old=c('CHR','FROM','TO'),new = range_keys)
          setkeyv(y,range_keys)
    }
)
names(all_snps) = c('denisova_file','neanderthal_file')

## separate naSNPs from aSNPs 
## and get the naSNPs that are in common between the two files

naSNPs = copy(all_snps)%>% lapply(function(x)x=x[ancestry=='non_archaic'][,c(11:14):= NULL])
naSNPs = naSNPs[[1]][naSNPs[[2]],on=c(range_keys,'REF',"ALT",'ANC','POP_ARCH_REF','POP_ARCH_ALT','POP_NOTARCH_REF','POP_NOTARCH_ALT','ancestry'),nomatch=0]
# nrow(naSNPs)
# [1] 5785079

## remove archaic/nonarchaic singletons
## remove those instances called as aSNPs in the other file
d_aSNPs = copy(all_snps[[1]])[!naSNPs,on=range_keys][ancestry != 'non_archaic']
# nrow(d_aSNPs)
# [1] 1339865
d_aSNPs <- d_aSNPs[,sum_asnp:=POP_ARCH_REF+POP_ARCH_ALT][sum_asnp>1][,sum_asnp:=NULL]
# nrow(d_aSNPs)
# 1] 1048497

n_aSNPs = copy(all_snps[[2]])[!naSNPs,on=range_keys][ancestry != 'non_archaic']
# nrow(n_aSNPs)
# [1] 966075
n_aSNPs <- n_aSNPs[,sum_asnp:=POP_ARCH_REF+POP_ARCH_ALT][sum_asnp>1][,sum_asnp:=NULL]
# nrow(n_aSNPs)
# [1] 756397

naSNPs <-naSNPs[,keep:=ifelse(ANC=='0' & POP_NOTARCH_ALT==1,'n',ifelse(ANC=='1' & POP_NOTARCH_REF==1,'no','yes'))][keep=='yes'][,keep:=NULL]
# nrow(naSNPs)
# [1] 4615313

## create single list with all sorted SNPs
## remove SNPs with unknown ANC 
## calculate MIAF/DAF and remove rare variants
## determine whether the main alleles are derived or ancestral
sorted_SNPs = list(d_aSNPs,n_aSNPs,naSNPs)%>%lapply(
  function(x){x=x[!ANC=='-1']%>%calculate_maf()
  }
)
names(sorted_SNPs) = c('denisova','neanderthal','modern_humans')
lapply(sorted_SNPs,function(x)nrow(x))
# $denisova
# [1] 1031260

# $neanderthal
# [1] 743025

# $modern_humans
# [1] 4476996

## plot SFS
sfs_before_filter <- copy(sorted_SNPs)%>%lapply(function(x){x=x[,MAF:=round(MAF,2)][,c(..range_keys,'MAF','ancestry')]})
sfs_before_filter <- Map(mutate,sfs_before_filter,ancestry=names(sfs_before_filter))%>%rbindlist()

## plot the results
pdf(paste(plot_dir,'sfs_before_filter.pdf',sep=''),width=12,height = 8)
ggplot(sfs_before_filter, aes(x=round(MAF,1),fill = ancestry))+ 
  geom_bar(colour="black")+
  scale_fill_manual(values = my_palette_pop)+
  xlab(' ')+ylab('Counts')+
  facet_wrap(~ancestry,ncol=3,scales='fixed')+
  theme_classic()+
  theme(
    legend.position='bottom'
  )
dev.off()

# sorted_snps_maf_fitered <- lapply(copy(sorted_SNPs),function(x){x=subset(x,MAF>=0.05)})

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
  ]
)
lapply(afr_filtered_snps,function(x)nrow(x))
# $denisova
# [1] 454277

# $neanderthal
# [1] 316412

# $modern_humans
# [1] 1108747

asnps_distr_btw_haps <- copy(afr_filtered_snps[c(1:2)])
asnps_distr_btw_haps <- Map(mutate,asnps_distr_btw_haps,ancestry=names(asnps_distr_btw_haps))%>%rbindlist()
asnps_distr_btw_haps <- asnps_distr_btw_haps[,snp_distribution_across_haps:=round(snp_distribution_across_haps,2)][,c('ancestry','snp_distribution_across_haps')]

## plot aSNP distribution across haplotypes results
pdf(paste(plot_dir,'asnps_distr_btw_haps.pdf',sep=''),width=12,height = 7)
ggplot(asnps_distr_btw_haps, aes(x=round(snp_distribution_across_haps,2),fill = ancestry))+ 
  geom_bar(colour="black")+
  scale_fill_manual(values = my_palette_pop)+
  xlab(' ')+ylab('Counts')+
  facet_wrap(~ancestry,ncol=3,scales='fixed')+
  theme_classic()+
  theme(
    legend.position='bottom'
  )
dev.off()

afr_filtered_snps <- lapply(afr_filtered_snps,function(x) x<-x[
  snp_distribution_across_haps >= 0.25
  ][
    ,snp_distribution_across_haps:=NULL
  ]
)
lapply(afr_filtered_snps,function(x)nrow(x))
# $denisova
# [1] 153336

# $neanderthal
# [1] 97399

# $modern_humans
# [1] 1108747

## plot SFS (2)
sfs_after_filter <- copy(afr_filtered_snps)%>%lapply(function(x){x=x[,MAF:=round(MAF,2)][,c(..range_keys,'MAF','ancestry')]})
sfs_after_filter <- Map(mutate,sfs_after_filter,ancestry=names(sfs_after_filter))%>%rbindlist()

## plot the results
pdf(paste(plot_dir,'sfs_after_filter.pdf',sep=''),width=12,height = 7)
ggplot(sfs_after_filter, aes(x=round(MAF,1),fill = ancestry))+ 
  geom_bar(colour="black")+
  scale_fill_manual(values = my_palette_pop)+
  xlab(' ')+ylab('Counts')+
  facet_wrap(~ancestry,ncol=3,scales='fixed')+
  theme_classic()+
  theme(
    legend.position='bottom'
  )
dev.off()

## remove ambiguous aSNPs, i.e. those in common between deni and nean
## try to assign snp to archaic ancestry based on freq within haps
## and then by allele state matching reference genome
## label all the remaining as ambiguous and remove them
ambiguous_asnps = afr_filtered_snps[[1]][
  afr_filtered_snps[[2]],on=c(range_keys,'REF',"ALT",'ANC','ancestry','DENI_REF','DENI_ALT','NEAN_REF','NEAN_ALT'),nomatch=0
  ][
    , assign := ifelse(round(MAF,2) > round(i.MAF,2), 'denisova',ifelse(round(MAF,2) == round(i.MAF,2),'ambiguous','neanderthal'))
]

ambiguous <- copy(ambiguous_asnps)[assign=='ambiguous'][
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

assigned_deni <- rbind(copy(ambiguous_asnps)[assign=='denisova'][,c(1:17)],copy(ambiguous)[match_genome=='denisova'][,c(1:17)])
assigned_nean <- rbind(copy(ambiguous_asnps)[assign=='neanderthal'][,c(1:17)],copy(ambiguous)[match_genome=='neanderthal'][,c(1:17)])
final_ambiguous <- copy(ambiguous)[match_genome=='ambiguous']

deni_asnps <- copy(afr_filtered_snps[[1]])[
  !rbind(copy(assigned_nean[,c(1:3)]),copy(final_ambiguous[,c(1:3)])),on=c(range_keys)
][,aoi:=ifelse(POP_ARCH_ALT>POP_ARCH_REF,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','aoi','state_allele')]

nean_asnps <- copy(afr_filtered_snps[[2]])[
  !rbind(copy(assigned_deni[,c(1:3)]),copy(final_ambiguous[,c(1:3)])),on=c(range_keys)
][,aoi:=ifelse(POP_ARCH_ALT>POP_ARCH_REF,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','aoi','state_allele')]

nasnps = copy(afr_filtered_snps[[3]])[,aoi:=ifelse(ANC==0,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','aoi','state_allele')]

filtered_snps = list(deni_asnps,nean_asnps,nasnps)
names(filtered_snps) = c('denisova','neanderthal','modern_humans')
lapply(filtered_snps,function(x)nrow(x))

# $denisova
# [1] 140916

# $neanderthal
# [1] 88625

# $modern_humans
# [1] 1108747

### combine asnps and subsample nasnps in order to match the archaic SFS
set.seed(2022)
counts_binned_asnps <- copy(rbind(filtered_snps$denisova,filtered_snps$neanderthal))[,roundMaf:=round(MAF,1)][,snp_bincount:=.N,by=.(roundMaf)][,c('roundMaf','snp_bincount')]%>%unique()
counts_binned_asnps <-rbind(counts_binned_asnps,data.table(roundMaf=1.0,snp_bincount=0))%>%setorderv('roundMaf',1)%>%split(by='roundMaf')

full_filtered_nasnps <- copy(filtered_snps$modern_humans)[,roundMaf:=round(MAF,1)]%>%setorderv('roundMaf',1)%>%split(by='roundMaf')

subsampled_nasnps <- purrr::map2(counts_binned_asnps,full_filtered_nasnps,function(x,y){
  subsampl<-sample_n(y, x$snp_bincount)
})%>%rbindlist()
subsampled_nasnps <- subsampled_nasnps[,roundMaf:=NULL]

final_set <- list(filtered_snps$denisova,filtered_snps$neanderthal,subsampled_nasnps)
final_set <- lapply(final_set,function(x)x<-x[,MAF:=round(MAF,1)][,allfreq:=ifelse(MAF<=0.05,'low',ifelse(MAF>0.2,'high','interm'))])
names(final_set) = c('denisova','neanderthal','modern_humans')

## plot the final sfs for archaics and non-archaics
pdf(paste(plot_dir,'sfs_archaics_nonarchaics.pdf',sep=''),width=12,height = 7)
df <-copy(final_set)
df <- Map(mutate,df,ancestry=names(df))%>%rbindlist()
df <- df[,ancestry:=ifelse(ancestry=='modern_humans','non_archaic','archaic')]
ggplot(df, aes(x=round(MAF,1),fill = ancestry))+ 
  geom_bar(colour="black")+
  scale_fill_manual(values = c('#81b29a','#f2cc8f'))+
  xlab(' ')+ylab('Counts')+
  facet_wrap(~ancestry,ncol=2,scales='fixed')+
  theme_classic()+
  theme(
    legend.position='bottom'
  )
dev.off()

## plot values for Bstatistic for the filtered set of SNPs
## do it for all snps and only those at highfreq
## for comparison add also all protein-coding sequences
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

all_exons <- exons(txdb)%>%as.data.table()
all_exons <- all_exons[!seqnames%like% "M|Un|_"][,c(..range_keys)][,ancestry:='all_exons']
setkeyv(all_exons,range_keys)

all_exons_bvals <- foverlaps(all_exons,bstatistic,type='any')[,c(..range_keys,'bvalue','ancestry')]%>%unique()%>%na.omit()

all_exons_col <- '#264653'
names(all_exons_col) = 'all_exons'

allsnps_bvalues <- copy(final_set)%>%lapply(function(x){
  overlaps <- foverlaps(x,bstatistic,type='within')[,c(..range_keys,'bvalue','MAF')]%>%unique()
})
allsnps_bvalues <- Map(mutate,allsnps_bvalues,ancestry=names(allsnps_bvalues))%>%rbindlist()

## plot the results
comparisons = list(
  c('denisova','modern_humans'),
  c('denisova','neanderthal'),
  c('neanderthal','modern_humans')
)

pdf(paste(plot_dir,'allsnps_bvalues.pdf',sep=''),width=12,height = 7)
snpsbval<- unique(copy(allsnps_bvalues)[,MAF:=NULL])
df<-rbind(snpsbval,all_exons_bvals)
ggplot(df, aes(x=ancestry,y=bvalue,fill =ancestry))+ 
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=F)+
  scale_fill_manual(values = c(my_palette_pop,all_exons_col))+
  xlab('Ancestry')+ylab('B-statistic')+
  stat_compare_means(
  method = "wilcox.test",
  comparisons=comparisons,
  size=5
  )+
  theme_classic()+
  theme(
    legend.position='none'
  )
dev.off()


pdf(paste(plot_dir,'highfreq_snps_bvalues.pdf',sep=''),width=12,height = 7)
snpsbval<- unique(copy(allsnps_bvalues)[MAF>0.05][,MAF:=NULL])
df<-rbind(snpsbval,all_exons_bvals)
ggplot(df, aes(x=ancestry,y=bvalue,fill =ancestry))+ 
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=F)+
  scale_fill_manual(values = c(my_palette_pop,all_exons_col))+
  xlab('Ancestry')+ylab('B-statistic')+
  stat_compare_means(
  method = "wilcox.test",
  comparisons=comparisons,
  size=5
  )+
  theme_classic()+
  theme(
    legend.position='none'
  )
dev.off()

## get gnomad allele frequencies
library(GenomicScores)
gnomadMAF <- MafDb.gnomAD.r3.0.GRCh38


filenames = paste0(out_dir,'new_',names(final_set),sep='.txt')
mapply(write.table,final_set, file = filenames,col.names = T, row.names = F, sep = "\t", quote = F)
