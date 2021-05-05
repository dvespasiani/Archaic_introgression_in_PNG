## Annotate SNPs
library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(annotatr)


setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')
input_dir = './filtered_files'
plot_dir='./Results/Plots/Annotation/'

options(width = 150)

range_keys = c('seqnames','start','end')

snps = list.files(input_dir,full.names = T,recursive = F) %>% lapply(function(y)fread(y,sep='\t',header = T))
names(snps)=c('denisova','modern_human','neanderthal')

## get genomic annotations
annots = c('hg19_basicgenes', 'hg19_genes_intergenic')
annotations = build_annotations(genome = 'hg19', annotations = annots)

## annotate SNPs
snps_annot = lapply(snps,function(y)
    y=makeGRangesFromDataFrame(y,keep.extra.columns=T)%>%
    annotate_regions(annotations = annotations,ignore.strand = T, quiet = F) %>% 
    as.data.table()
    )%>%
    lapply(function(y)y=y[
            ,genomic_element:=plyr::revalue(
                annot.type,c(
                    'hg19_genes_exons'='Exons',
                    'hg19_genes_introns'='Introns',
                    'hg19_genes_intergenic'='Intergenic',
                    'hg19_genes_promoters'='CRE',
                    'hg19_genes_1to5kb'='CRE',
                    'hg19_genes_3UTRs'='UTRs',
                    'hg19_genes_5UTRs'='UTRs')
                    )
                    ][
                        ,c(..range_keys,'REF',"ALT",'MAF','genomic_element')
                        ]%>%unique()
)

## calculate odds ratio SNPs in genomic elements 
counts_snps_annot = copy(snps_annot)%>%
    lapply(function(x)x=x[
            ,numb_snps:= .N,
            ][
                ,numbsnps_element:= .N, by=.(genomic_element)
                ][
                    ,numbsnps_non_element := numb_snps-numbsnps_element 
                    ][
                        ,c('genomic_element','numbsnps_element','numbsnps_non_element')
                        ]%>%unique()
)
counts_snps_annot = Map(mutate,counts_snps_annot,file = names(counts_snps_annot))

calculate_odds_ratio = function(asnps,nasnps){
    table = asnps[
        nasnps,on='genomic_element',nomatch=0
        ][
            ,c('file','i.file'):=NULL
        ]%>%
        split(by='genomic_element')%>%
        lapply(function(x)x=x[,genomic_element:=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%fisher.test())
        
    fisher_test_results = copy(table)%>%
    lapply(function(x)x=data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
    ))
    fisher_test_results = Map(mutate,fisher_test_results,elements=names(fisher_test_results))%>%rbindlist()
    fisher_test_results = fisher_test_results[,significance:=ifelse(p<0.05,'*','')]

    return(fisher_test_results)
}

deni_odds_ratio = calculate_odds_ratio(counts_snps_annot[[1]],counts_snps_annot[[2]])[,ancestry:='denisova']
nean_odds_ratio = calculate_odds_ratio(counts_snps_annot[[3]],counts_snps_annot[[2]])[,ancestry:='neanderthal']

asnps_odds_ratio = rbind(deni_odds_ratio,nean_odds_ratio)

## plot the results
palette =c(
    'thistle3', #intergenic
    'olivedrab3', #intron
    'lightgoldenrod', # CRE
    'tan3', #exons
    'darkorchid1' # utr
)
names(genomic_element_palette) = unique(asnps_odds_ratio$elements)

pdf(paste(plot_dir,'odds_ratio_asnps_vs_nasnps.pdf',sep=''),width=8,height = 5)
ggplot(asnps_odds_ratio, aes(x=elements, y=odds_ratio,label = significance)) + 
geom_point(aes(colour = elements))+
geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
scale_colour_manual(values = genomic_element_palette)+
geom_hline(yintercept=1,linetype='dashed',size=.5)+
geom_text(aes(y= asnps_odds_ratio$upper_ci+0.1),size=5)+
xlab(' ')+ylab('odds ratio \n aSNPs vs naSNPs')+
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
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()





