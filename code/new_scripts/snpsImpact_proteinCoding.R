
library(data.table)
library(magrittr)
library(dplyr)
library(GenomicRanges)
library(annotatr)
library(R.utils)
library(ggthemes)
library(ggplot2)
library(ggpubr)


options(width = 150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

input_dir = './filtered_files'
plot_dir = './Results/Plots/genomeAnnot/'

snps = list.files(input_dir,full.names = T,recursive = F,pattern='new_*') %>% lapply(function(y)fread(y,sep='\t',header = T))
names(snps)=ancestry

## annote variants to find those within protein coding sequences
annots <- c('hg19_basicgenes')

annotations = build_annotations(genome = 'hg19', annotations = annots)

snpsAnnot <- copy(snps)%>%lapply(function(x){
    annoX = annotate_regions(
    regions = makeGRangesFromDataFrame(x,keep.extra.columns=T),
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)%>%as.data.table()
    annoX <- annoX%>%dplyr::select(all_of(names(x)),'annot.type')%>%unique()
    return(annoX)
})

## probability of aSNPs lying within protein-coding sequences 
## is calculated as the tot numb of human protein coding snps/tot numb human snps

numbProteinCodingSNPs <- copy(snpsAnnot[[2]])[annot.type=='hg19_genes_exons'][,c(..range_keys)]%>%unique()%>%nrow()
numbSNPs <- copy(snpsAnnot[[2]])[,c(..range_keys)]%>%unique()%>%nrow()
probProteinCoding <- numbProteinCodingSNPs/numbSNPs

numbSuccesses<-copy(snpsAnnot[c(1,3)])%>%lapply(function(x)x[annot.type=='hg19_genes_exons'][,c(..range_keys)]%>%unique()%>%nrow())
numbTrials<-copy(snpsAnnot[c(1,3)])%>%lapply(function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

binPval <-purrr::map2(numbSuccesses,numbTrials,function(x,n){
    binomialRes <-data.table(
        pval =binom.test(x,n,probProteinCoding)$p.value,
        numbSuccess = binom.test(x,n,probProteinCoding)$statistic
    )
    return(binomialRes)
})
# $denisova
#            pval numbSuccess
# 1: 1.413473e-11        3431

# $neanderthal
#            pval numbSuccess
# 1: 6.637236e-19        1919

## proportion variants per genomic elements
allSNPsAnnot <- copy(snpsAnnot)%>%lapply(function(x){
    x<-x[,numbSNPsAnnot:=.N,by=.(annot.type)][,numbSNPs:=.N][,propSNPsAnnot:=round((numbSNPsAnnot/numbSNPs)*100,2)]
})

commonToHighFreq_SNPsAnnot <- copy(snpsAnnot)%>%lapply(function(x){
    x<-x[allfreq!='low'][,numbSNPsAnnot:=.N,by=.(annot.type)][,numbSNPs:=.N][,propSNPsAnnot:=round((numbSNPsAnnot/numbSNPs)*100,2)]
})

plot_propSNPsGenome <- function(x){
    df <- copy(x)
    df <- Map(mutate,df,ancestry=names(df))%>%rbindlist()
    df <- df[,c('ancestry','propSNPsAnnot','annot.type')]%>%unique()
    df <- df[,annot.type:=stringr::str_remove(annot.type, 'hg19_genes_')]
    p <- ggplot(df,aes(x=reorder(annot.type,-propSNPsAnnot),y=propSNPsAnnot,fill=ancestry))+
        geom_bar(stat='identity',position='dodge')+
        ylab('Proportion SNPs') + xlab('Genomic feature')+
        scale_fill_manual(values=my_palette_pop)+
        theme_classic()+
        theme(
            legend.position = 'bottom'
        )
    return(p)
}

pdf(paste(plot_dir,'GenomicFeature_allSNPs.pdf',sep=''),width=7,height = 7)
plot_propSNPsGenome(allSNPsAnnot)
dev.off()


pdf(paste(plot_dir,'GenomicFeature_commonTohighfreq.pdf',sep=''),width=7,height = 7)
plot_propSNPsGenome(commonToHighFreq_SNPsAnnot)
dev.off()