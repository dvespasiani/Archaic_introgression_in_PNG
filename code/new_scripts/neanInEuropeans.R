## script used to control whether pipeline works also for neanderthal snps in europeans
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(GenomicRanges)

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
neanSnpsInEur <- fread(list.files(input_dir,recursive =F,full.names = T,pattern='europe'),sep='\t')[
        POP_ARCH_REF + POP_ARCH_ALT>0,ancestry := 'archaic' 
        ][
            is.na(ancestry), ancestry := 'non_archaic'
]%>%setnames(old=c('CHR','FROM','TO'),new = range_keys)
setkeyv(neanSnpsInEur,range_keys)


## separate naSNPs from aSNPs 
## and get the naSNPs that are in common between the two files
naSNPs <- copy(neanSnpsInEur)[ancestry=='non_archaic'][,c(11:14):= NULL]

numbSnps <- nrwo(neanSnpsInEur)
# [1] 11422574
numbNASnps <- nrwo(naSNPs)
# [1] 9356675

aSNPs <- copy(neanSnpsInEur)[!naSNPs,on=range_keys][ancestry != 'non_archaic'][,c(11:14):= NULL]
# nrow(aSNPs)
# [1] 2065899

aSNPsNoSingletons <- copy(aSNPs)[,sum_asnp:=POP_ARCH_REF+POP_ARCH_ALT][sum_asnp>1][,sum_asnp:=NULL]
# nrow(aSNPsNoSingletons)
# [1] 1496642

naSNPsNoSingletons <- copy(naSNPs)[,keep:=ifelse(ANC=='0' & POP_NOTARCH_ALT==1,'n',ifelse(ANC=='1' & POP_NOTARCH_REF==1,'no','yes'))][keep=='yes'][,keep:=NULL]
# nrow(naSNPsFiltered)
# [1] 6421584

## create single list with all sorted SNPs
## remove SNPs with unknown ANC 
## calculate MIAF/DAF and remove rare variants
## determine whether the main alleles are derived or ancestral
sortedSNPs = list(aSNPsNoSingletons,naSNPsNoSingletons)%>%lapply(
  function(x){x=x[!ANC=='-1']%>%calculate_maf()
  }
)
names(sortedSNPs) = c('neanderthal','modern_humans')
lapply(sortedSNPs,function(x)nrow(x))
# $neanderthal
# [1] 1470213

# $modern_humans
# [1] 6211120

## remove SNPs segregating within africans
afr_1kg_snp_freq  = fread('../Annotation_and_other_files/human_genome/1kg_snp_modified.gz',sep=' ',header=T,select=c(1:5))
colnames(afr_1kg_snp_freq)[1:2] = range_keys[-3]

snps_notin_afr = lapply(sortedSNPs,function(x) x=x[afr_1kg_snp_freq,on=c(range_keys[-3],'REF','ALT'),nomatch=0][AFR_AF<=0.005][,AFR_AF:=NULL])
snps_notin_1kg = lapply(sortedSNPs,function(x)x=x[!afr_1kg_snp_freq,on=c(range_keys[-3],'REF','ALT')])

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

asnps_distr_btw_haps <- copy(afr_filtered_snps)
asnps_distr_btw_haps <- Map(mutate,asnps_distr_btw_haps,ancestry=names(asnps_distr_btw_haps))%>%rbindlist()
asnps_distr_btw_haps <- asnps_distr_btw_haps[,snp_distribution_across_haps:=round(snp_distribution_across_haps,2)][,c('ancestry','snp_distribution_across_haps')]


afr_filtered_snps <- lapply(afr_filtered_snps,function(x) x<-x[
  snp_distribution_across_haps >= 0.25
  ][
    ,snp_distribution_across_haps:=NULL
  ]
)

neanASNPs <- copy(afr_filtered_snps[[1]])[
    ,aoi:=ifelse(POP_ARCH_ALT>POP_ARCH_REF,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','aoi','state_allele')
]

humanNASNPs = copy(afr_filtered_snps[[2]])[
    ,aoi:=ifelse(ANC==0,'alt','ref')][,c(..range_keys,"REF",'ALT','MAF','aoi','state_allele')
]

filtered_snps = list(neanASNPs,humanNASNPs)
names(filtered_snps) = c('neanderthal','modern_humans')
lapply(filtered_snps,function(x)nrow(x))
# $neanderthal
# [1] 118347

# $modern_humans
# [1] 1502154

### combine asnps and subsample nasnps in order to match the archaic SFS
set.seed(2022)
counts_binned_asnps <- copy(filtered_snps$neanderthal)[,roundMaf:=round(MAF,1)][,snp_bincount:=.N,by=.(roundMaf)][,c('roundMaf','snp_bincount')]%>%unique()
counts_binned_asnps <-rbind(counts_binned_asnps,data.table(roundMaf=1.0,snp_bincount=0))%>%setorderv('roundMaf',1)%>%split(by='roundMaf')

full_filtered_nasnps <- copy(filtered_snps$modern_humans)[,roundMaf:=round(MAF,1)]%>%setorderv('roundMaf',1)%>%split(by='roundMaf')
full_filtered_nasnps <- full_filtered_nasnps[names(full_filtered_nasnps) %in% names(counts_binned_asnps)]

subsampled_nasnps <- purrr::map2(counts_binned_asnps,full_filtered_nasnps,function(x,y){
  subsampl <- sample_n(y, x$snp_bincount)
})%>%rbindlist()
subsampled_nasnps <- subsampled_nasnps[,roundMaf:=NULL]

final_set <- list(filtered_snps$neanderthal,subsampled_nasnps)
final_set <- lapply(final_set,function(x)x<-x[,MAF:=round(MAF,1)][,allfreq:=ifelse(MAF<=0.05,'low',ifelse(MAF>0.2,'high','interm'))])
names(final_set) = c('neanderthal','modern_humans')

##===============================
## RECALL NEAN SNPS IN EUR 
##===============================

# check how many SNPs from danneman et al 2017 you recover
# knownNeanSNPs <- list.files('./knownNeanInEUR/',recursive=F,full.names=T,pattern='neanSNPseQTLs_findley2021.txt')%>%
# lapply(function(x){
#     snps <- fread(x,sep='\t',header=T)[
#         ,seqnames:=paste('chr',chr,sep='')
#         ][
#             ,chr:=NULL
#             ][
#                 ,start:=pos
#                         ][
#                             ,pos:=NULL
#                             ]
#     return(snps)
# })
# names(knownNeanSNPs) = gsub(".*_","",gsub('\\.txt*','',list.files('./knownNeanSNPs/',recursive=F,full.names=F)))


# neanGO = great_enrichment(cresSnps[[1]])

# myNeanGenes <- copy(neanGO[[2]])[,c('gene')]%>%unique()%>%na.omit()

# ## list of knonw genes targeted by nean snps in eur
# knownNeanGenes <- list.files('./knownNeanInEUR',recursive=F,full.names=T,pattern='*Genes*')%>%
# lapply(function(x){ x<- fread(x,sep='\t',header=T,stringsAsFactors = F)%>%na.omit()})

# names(knownNeanGenes) = gsub(".*_","",gsub('\\.txt*','',list.files('./knownNeanInEUR',recursive=F,full.names=F,pattern='*Genes*')))

# ## make upset plot
# library(ComplexHeatmap)

# listNeanGenes <- c(knownNeanGenes,list(myNeanGenes))%>%lapply(function(x)(x$gene))
# names(listNeanGenes) = make.names(c(names(knownNeanGenes),'myNeanGenes'))

# upsetMatrix = make_comb_mat(list_to_matrix(listNeanGenes))

# colors = c('#e76f51','#2a9d8f','#e9c46a')
# names(colors) = names(listNeanGenes)

# pdf(paste(plot_dir,'upsetMyRecalledNeanGenes.pdf',sep=''),width=7,height = 5)
# UpSet(
#     upsetMatrix,
#     right_annotation = upset_right_annotation(
#         upsetMatrix, 
#         gp = gpar(fill = colors)
#     )
# )
# dev.off()



# commonNeanGenes <- data.table(gene=Reduce(intersect, listNeanGenes))



# dannemann2017NeanSNPs <- copy(knownNeanSNPs[[4]])[
#     ,ALT:=stringr::str_replace_all(ALT, "[[:punct:]]", "")
#     ][
#         ,REF:=stringr::str_replace_all(REF, "[[:punct:]]", "")
#         ][
#             ,MAF:=as.numeric(neanAlleleFreq)
#         ][
#             ,neanAlleleFreq:=NULL
# ]
# silvert2019NeanSNPs <- copy(knownNeanSNPs[c(1,3)])%>%lapply(function(x){
#     x<- x[,MAFCEU:=as.numeric(stringr::str_replace_all(MAFCEU, "[[:punct:]]", ""))][,MAF:=MAFCEU/1000][,MAFCEU:=NULL]
# })%>%rbindlist()%>%unique()

# findley2021NeanSNPs <- copy(knownNeanSNPs[[2]])[,MAF:=Nean_AF]

# finalKnownNeanSNPs <- list(silvert2019NeanSNPs,findley2021NeanSNPs,dannemann2017NeanSNPs)
# names(finalKnownNeanSNPs) = names(knownNeanSNPs)[-1]

# numbKnownNeanSNP <- lapply(finalKnownNeanSNPs,function(x)nrow(x[,c(..range_keys[-3])]%>%unique()))

# recallNeanNoFilter <- copy(finalKnownNeanSNPs)%>%lapply(function(x){
#     x<-x[aSNPsNoSingletons,on=c(range_keys[-3]),nomatch=0]
# })

# numbRecallNeanSNPNoFilt <- lapply(recallNeanNoFilter,function(x)nrow(x[,c(..range_keys)]%>%unique()))

# filteredSetNeanSNPs <- copy(final_set[[1]])

# recallNeanFiltered <- copy(finalKnownNeanSNPs)%>%lapply(function(x){
#     x<-x[filteredSetNeanSNPs,on=c(range_keys[-3]),nomatch=0]
# })

# numbRecallNeanSNPFilt <- lapply(recallNeanFiltered,function(x)nrow(x[,c(..range_keys)]%>%unique()))

# propRecallNoFilt <- purrr::map2(numbRecallNeanSNPNoFilt,numbKnownNeanSNP,function(x,y){round((x/y)*100,2)})

# propRecallFilt <- purrr::map2(numbRecallNeanSNPFilt,numbKnownNeanSNP,function(x,y){round((x/y)*100,2)})

# testCorr <- lapply(recallNeanFiltered,function(x)cor.test(x$MAF,x$i.MAF))

# ##===================
# ## make upset plot
# ##===================
# library(ComplexHeatmap)

# getSNPsIDs <- function(x){
#     df<- copy(x)%>%lapply(function(y){
#         y<-y[,snpID:=paste(seqnames,start,sep='_')][,snpID]%>%unique()
#         })
    
    
#     return(df)
# }

# recalledSNPsIDs <- getSNPsIDs(recallNeanFiltered)

# upsetMatrix = make_comb_mat(list_to_matrix(recalledSNPsIDs))

# colors = c('#e76f51','#2a9d8f','#e9c46a')
# names(colors) = names(recallNeanFiltered)

# pdf(paste(plot_dir,'upsetMyRecalledNeanSNPs.pdf',sep=''),width=7,height = 5)
# UpSet(
#     upsetMatrix,
#     right_annotation = upset_right_annotation(
#         upsetMatrix, 
#         gp = gpar(fill = colors)
#     )
# )
# dev.off()

# ##==============================
# ## make proportion recall plot 
# ##==============================
# makeDT <- function(list,type){
#     dt <-copy(list)%>%lapply(function(x){x<-data.table(prop=x,type=type)})
#     dt <- Map(mutate,dt,file=names(dt))%>%rbindlist()
#     return(dt)
# }
# propRecallFiltDT <- makeDT(propRecallFilt,'filtered')
# propRecallNoFiltDT <- makeDT(propRecallNoFilt,'no_filtered')

# propDT <- rbind(propRecallFiltDT,propRecallNoFiltDT)

# pdf(paste(plot_dir,'propMyRecalledNeanSNPs.pdf',sep=''),width=7,height = 5)
# ggplot(propDT,aes(x=file,y=prop,fill=file))+
# geom_bar(stat='identity')+
# xlab(' ')+ylab('SNP recall proportion')+
# scale_fill_manual(values=colors)+
# facet_wrap(type~.,ncol=2)+
# theme_classic()+
# theme(
#     legend.position='bottom'
# )
# dev.off()

# ##==============================
# ## make MAF correlation plot
# ##==============================
# recallNeanFiltered_corrMAF <- copy(recallNeanFiltered)%>%lapply(function(x){
#     x<-x[,snpID:=paste(seqnames,start,sep='_')][,c('MAF','i.MAF','snpID')]%>%setnames(c('myMAF','reportedMAF','snpID'))
# })
# recallNeanFiltered_corrMAF <- Map(mutate,recallNeanFiltered_corrMAF,file=names(recallNeanFiltered_corrMAF))%>%rbindlist()

# # mafCorrMatrix <- copy(recallNeanFiltered_corrMAF)
# # rowNames = mafCorrMatrix$snpID


# # pdf(paste(plot_dir,'mafCorrelationMyRecalledNeanSNPs.pdf',sep=''),width=7,height = 5)
# # ggscatter(
# #     recallNeanFiltered_corrMAF, x = "myMAF", y = "reportedMAF",
# #           add = "reg.line",                                 
# #           conf.int = TRUE,  
# #           facet.by = 'file',  
# #           color = "file", palette = colors,                   
# #           add.params = list(color = "blue",fill = "lightgray")
# #           )+
# #           stat_cor(method = "pearson")+
# #           theme_classic()
# # dev.off()

# vector_of_scores <- c(testCorr[[1]]$estimate, testCorr[[2]]$estimate,testCorr[[3]]$estimate)
# matrix_of_scores <- matrix(vector_of_scores,3,1)
# rownames(matrix_of_scores) = names(testCorr)
# colnames(matrix_of_scores) = 'my'
# melted_matrix = reshape2::melt(matrix_of_scores,na.rm=T)

# pdf(paste(plot_dir,'mafCorrelationMyRecalledNeanSNPs.pdf',sep=''),width=10,height = 7)
# p <- ggplot(melted_matrix, aes(Var2, Var1, fill = value))+
#  geom_tile(color = "white")+
#  scale_fill_gradient2(
#    low=c("blue"),mid=c("white"),high=c("red"),
#    midpoint = 0, limit = c(-1,1), space = "Lab", 
#    name="Pearson correaltion") +
#   theme_minimal()+ # minimal theme
#   theme(
#     axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)
#     )+
#     coord_fixed()
# p + geom_text(
#   aes(Var2, Var1, label = value), 
#   color = "black", size = 4) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = 'bottom',
#     legend.direction = "horizontal"
#     )+
#     guides(
#       fill = guide_colorbar(barwidth = 7, barheight = 1,
#       title.position = "top", title.hjust = 0.5)
# )
# dev.off()



## now annotate with chromHMM and look enrichment over chromatin states and tissues
cells_dir = '../Annotation_and_other_files/Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states'

## load cells
cells = list.files(cells_dir,recursive = F,full.names = T) %>% 
    lapply(function(y)fread(y,sep = ' ',header=T) %>% makeGRangesFromDataFrame(keep.extra.columns =T) %>% as.data.table()
)
cells = lapply(cells,function(x)x=x[,width:=as.numeric(width)][,state_coverage:=sum(width),by=.(chrom_state,cell_type)])
names(cells) = gsub("\\..*","",list.files(cells_dir,recursive=F,full.names=F))
cells = cells[c(1:18)] # remove the encode ones

### merge snps with chromatin states
lapply(final_set,function(x)setkeyv(x,range_keys))
lapply(cells,function(x)setkeyv(x,range_keys))

chromatin_state_annotation = function(cell,alleles){
  cell=lapply(cell,function(cell)
    alleles=lapply(alleles,function(alleles)
      foverlaps(cell,alleles, type="within")[
        ,element_seqnames:=seqnames
      ][
        ,c('width','strand'):=NULL] %>%
        setnames(
          old= c(range_keys[-1],paste('i.',range_keys[-1],sep='')),
          new = c(paste('element',range_keys[-1],sep='_'),range_keys[-1])
          )
    )
  )
}

snps_chromatin_states = chromatin_state_annotation(final_set,cells) 

neandertal_chromatin_states = rbindlist(snps_chromatin_states[[1]])
human_chromatin_states = rbindlist(snps_chromatin_states[[2]])

snps_chromatin_states = list(neandertal_chromatin_states,human_chromatin_states) %>% lapply(function(x)group_celltypes(x))
names(snps_chromatin_states) = c('neanderthal','modern_humans')

count_snps_chromstate <- function(snps,select_cols){
  if('cell_line' %in% select_cols){
  counts <- copy(snps)
  counts<-counts[
      ,numbsnps_chromstate:=.N,by=.(chrom_state,cell_type)
      ][
        ,numbsnps:=.N,by=.(cell_type)
        ][
          ,numbsnps_notchromstate:=numbsnps-numbsnps_chromstate
          ]
          
  } else{
  counts <- copy(snps)
  counts<-counts[
      ,numbsnps_chromstate:=.N,by=.(chrom_state)
      ][
        ,numbsnps:=.N,
        ][
          ,numbsnps_notchromstate:=numbsnps-numbsnps_chromstate
          ]
         
  }
  final_counts <- copy(counts)%>%dplyr::select(c('numbsnps_chromstate','numbsnps_notchromstate',all_of(select_cols)))%>%unique()
  return(final_counts)
}

split_element_col <- function(x){
  x=copy(x)[,c('chrom_state','cell_type','cell_line'):= tstrsplit(elements, ".", fixed=TRUE)]
  return(x)
}

## genome wide enrichment (no cell type resolution)
snp_annot_counts <- copy(snps_chromatin_states)%>%lapply(function(x){
  x<-copy(x)[,c(..range_keys,'chrom_state')]%>%unique()%>%count_snps_chromstate(c('chrom_state'))
})

genwideNeanOR = calculate_or(snp_annot_counts[[1]],snp_annot_counts[[2]],c('chrom_state'))[,ancestry:='neanderthal']%>%adjust_pvalues()

plot_or <- function(or){
  p<-ggplot(or, aes(x=factor(elements, level = chrom_states), y=odds_ratio,label = p_signif)) + 
    geom_point(aes(colour = elements))+
    geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
    scale_colour_manual(values = chromHMM_palette)+
    geom_hline(yintercept=1,linetype='dashed',size=.5)+
    geom_text(aes(y= or$upper_ci+0.1),size=5)+
    xlab(' ')+ylab('OR aSNPs vs naSNPs')+
    theme_classic()+
    theme(
            panel.background =element_rect(fill = 'white', colour = 'black',size=1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.y = element_text(hjust = 0.5),
            strip.background = element_rect(color = 'black', linetype = 'solid'),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
            legend.position = "bottom",
            legend.key = element_rect(fill = "white", colour = "black"),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
            )
  return(p)
}

pdf(paste(plot_dir,'genwideASNPsOR_NeanInEUR.pdf',sep=''),width=5,height = 5)
plot_or(genwideNeanOR)
dev.off()

meansnp_annot_counts <- copy(snps_chromatin_states[[1]])[
    ,c(..range_keys,'chrom_state','cell_line','cell_type')]%>%
    unique()%>%
    count_snps_chromstate(c('cell_line','cell_type','chrom_state'))
meansnp_annot_counts<- meansnp_annot_counts[,numbsnps_notchromstate:=NULL][,meanNumbSnps:=round(mean(numbsnps_chromstate),1),by=.(chrom_state)
  ][
    ,c('chrom_state','meanNumbSnps')
]%>%unique()


pdf(paste(plot_dir,'genwideNumbASNPs_NeanInEUR.pdf',sep=''),width=10,height = 3)
ggplot(meansnp_annot_counts, aes(x=factor(chrom_state, level = chrom_states), y=log10(meanNumbSnps),fill=chrom_state)) + 
   geom_bar(stat='identity')+
    scale_fill_manual(values = chromHMM_palette)+
    xlab('')+ylab('log10 mean number aSNPs')+
    theme_classic()+
    theme(
            panel.background =element_rect(fill = 'white', colour = 'black',size=1),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.y = element_text(hjust = 0.5),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
            legend.position = "none",
            legend.key = element_blank(),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
            )
dev.off()



## repeat it with common-to-highfreq in cres
cresSnps <- copy(snps_chromatin_states)%>%lapply(function(x)x[chrom_state%in% cres_states][allfreq!='low'][,c(..range_keys,'cell_line')]%>%unique())

cresSnps = lapply(
  cresSnps,function(x)
  x=x[
    ,numbsnps:= .N
    ][
      ,numbsnps_in_tissue := .N,by=.(cell_line)
      ][
        ,numbsnps_notin_tissue := numbsnps - numbsnps_in_tissue
        ][
          ,c('cell_line','numbsnps_in_tissue','numbsnps_notin_tissue')
          ]%>%unique()
)

tissueNeanOR = calculate_or(cresSnps[[1]],cresSnps[[2]],c('cell_line'))[,ancestry:='neanderthal']%>%adjust_pvalues()

pdf(paste(plot_dir,'enrichmentCREsCommonToHighFreqASNPsOR_NeanInEUR.pdf',sep=''),width=5,height = 5)
ggplot(tissueNeanOR, aes(x=elements, y = odds_ratio,label = p_signif))+ 
  geom_point(aes(colour = elements))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
  scale_colour_manual(values = tissue_colors)+
  geom_hline(yintercept=1,linetype='dashed',size=.5)+
  geom_text(aes(y= tissueNeanOR$upper_ci+0.05),size=5)+
  xlab(' ')+ylab('OR aSNPs vs naSNPs')+
  facet_wrap(~ancestry)+
  theme(
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
          )
dev.off()


meanNumbaSNPS <- copy(snps_chromatin_states[[1]])[
      chrom_state %in% cres_states & allfreq!='low' 
      ][
    ,c(..range_keys,'cell_line','cell_type')]%>%
unique()
meanNumbaSNPS <- meanNumbaSNPS[
      ,numbsnps_celltype:=.N,by=.(cell_type)
      ][,meanNumbSnps:=round(mean(numbsnps_celltype),1),by=.(cell_line)
  ][
    ,c('cell_line','meanNumbSnps')
]%>%unique()


pdf(paste(plot_dir,'enrichmentCREsCommonToHighFreqNumbASNPs_NeanInEUR.pdf',sep=''),width=8,height = 3)
ggplot(meanNumbaSNPS, aes(x=tissue, y=log10(meanNumbSnps),fill=tissue)) + 
  geom_bar(stat='identity')+
  scale_fill_manual(values = tissue_colors)+
  xlab('')+ylab('log10 mean number aSNPs')+
  theme_classic()+
  theme(
        panel.background =element_rect(fill = 'white', colour = 'black',size=1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(hjust = 0.5),
        strip.background.y = element_blank(),
        strip.background.x =element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
)
dev.off()



## KLF3 

# ## plot numb snps per chrom state
# numbasnps_chromstate <- copy(snp_annot_counts[[1]])

# plot_numbsnps <- function(x){
#    p<-ggplot(x, aes(x=factor(chrom_state, level = chrom_states), y=log10(numbsnps_chromstate),fill=chrom_state)) + 
#    geom_bar(stat='identity')+
#     scale_fill_manual(values = chromHMM_palette)+
#     xlab('')+ylab('log10 number aSNPs')+
#     theme_classic()+
#     theme(
#             panel.background =element_rect(fill = 'white', colour = 'black',size=1),
#             panel.grid.minor = element_blank(),
#             panel.border = element_blank(),
#             panel.grid.major = element_blank(),
#             strip.text.y = element_text(hjust = 0.5),
#             strip.background.y = element_blank(),
#             strip.background.x =element_blank(),
#             legend.position = "none",
#             legend.key = element_blank(),
#             axis.line = element_blank(),
#             axis.ticks.x = element_blank(),
#             axis.text.x = element_blank()
#             )
#   return(p)
# }

# pdf(paste(plot_dir,'genwideNumbASNPs_NeanInEUR.pdf',sep=''),width=5,height = 3)
# plot_numbsnps(numbasnps_chromstate)
# dev.off()



