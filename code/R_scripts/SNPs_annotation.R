## Annotate SNPs
library(dplyr); library(data.table)
library(magrittr); 
library(purrr)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(annotatr)


numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

plot_dir='./Results/Plots/Annotation/'

snps=as.character(list.files('./Grouped_filtered_snps/new_set/',full.names = T,recursive = F)) %>% 
  lapply(function(y)fread(y,sep=' ',header = T) %>% 
           setnames(old=c('CHR','FROM','TO'),new=c('seqnames','start','end')))

names(snps)=c('denisova','neandertal','png')

annotation_function = function(x){
  annots = c('hg19_basicgenes', 'hg19_genes_intergenic')
  annotations = build_annotations(genome = 'hg19', annotations = annots)
  
  df=copy(x)
  df=lapply(df,function(y)
    y=y[,c('seqnames','start','end','all_frequency')] %>% unique() %>%
      makeGRangesFromDataFrame(keep.extra.columns=T))
  
  df=lapply(df,function(y)  
    y=annotate_regions(regions =y,
                       annotations = annotations,
                       ignore.strand = T, quiet = F) %>% as.data.table())
  df=lapply(df,function(y)y=y[,c('seqnames','start','end','annot.type','all_frequency')] %>% unique())
  return(df)
}

snps_annotated=annotation_function(snps)

##------------------------------------------------------------------------------------------------
## calculate the proportion of SNPs annotated for each genomic element per allele frequency bin
##------------------------------------------------------------------------------------------------

snps_annotated_proportion=copy(snps_annotated)

snps_annotated_proportion=lapply(snps_annotated_proportion,function(y)y=y[
  ,all_frequency:=round(all_frequency,1)
  ][
    ,genomic_element:=plyr::revalue(annot.type,c('hg19_genes_exons'='Exons',
                                                 'hg19_genes_introns'='Introns',
                                                 'hg19_genes_intergenic'='Intergenic',
                                                 'hg19_genes_promoters'='Promoters',
                                                 'hg19_genes_1to5kb'='Enhancers',
                                                 'hg19_genes_3UTRs'='UTRs',
                                                 'hg19_genes_5UTRs'='UTRs'))
    ][
      ,annot.type:=NULL
      ][
    ,numb_snps_freqbin:= .N, by=all_frequency
    ][
      ,numbsnps_element_freqbin:= .N, by=.(all_frequency,genomic_element)
      ][
        ,proportion_per_freqbin:= numbsnps_element_freqbin/numb_snps_freqbin
        ][
          ,log10numbsnps_per_freqbin:=log10(numb_snps_freqbin) 
          ])

snps_annotated_proportion=Map(mutate,snps_annotated_proportion,'pop'=names(snps_annotated_proportion))
snps_annotated_proportion=lapply(snps_annotated_proportion,function(x)x=setDT(x))
## plot proportion snps per element and freq bin + number of snps per freq bin

create_plots=function(x,pop_palette){
  proportion_element=copy(x)[,c('all_frequency','genomic_element','pop','proportion_per_freqbin')] %>% unique()
  number_snps=copy(x)[,c('log10numbsnps_per_freqbin','all_frequency','pop')] %>% unique()

  barcol=scale_fill_manual(name= " ",values = pop_palette)
  
  
  my_palette_genomicelement=c('lightgoldenrod', # Enhancers
                              'tan3', #exons
                              'thistle3', #intergenic
                              'olivedrab3', #intron
                              'orangered2', # Promoters
                              'darkorchid1' # utr
  )
  
  names(my_palette_genomicelement)=levels(as.factor(df$genomic_element))
  colScale_genomicelement=scale_fill_manual(name= " ",values = my_palette_genomicelement,
                                            labels = c('Enhancers',"Exon",
                                                       "Intergenic","Intron",
                                                       "Promoters",'UTR'))
  
  
  
 p1=ggplot(number_snps,aes(x=all_frequency,y=log10numbsnps_per_freqbin,fill=pop))+
   geom_bar(stat='identity',position=position_dodge())+
   barcol+
   xlab('')+
   ylab('\n Log10 number of SNPs \n')+
   theme(legend.position = 'none',
         panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line = element_blank())
 
 
  p2=ggplot(proportion_element,aes(x=all_frequency,y=proportion_per_freqbin,fill=genomic_element))+
    geom_bar(stat = 'identity',col='black')+
    colScale_genomicelement+
    ylab('\n Proportion of SNPs annotated \n')+ xlab('\n Allele frequency in PNG \n')+
    theme(legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"),
          panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_blank())
 
   p3=ggpubr::ggarrange(p1,p2,ncol=1,align = 'v')  
  return(p3)
}

pdf(paste(plot_dir,'Denisova_annotation.pdf',sep=''),width = 7,height = 7)
create_plots(snps_annotated_proportion[[1]],'#C99E10')
dev.off()

pdf(paste(plot_dir,'Neanderthal_annotation.pdf',sep=''),width = 7,height = 7)
create_plots(snps_annotated_proportion[[2]],'#9B4F0F')
dev.off()

pdf(paste(plot_dir,'PNG_annotation.pdf',sep=''),width = 7,height = 7)
create_plots(snps_annotated_proportion[[3]],'#1E656D')
dev.off()



## check if there is no statistical enrichment/depletion within exons

exons=function(x){
  df=copy(x)
  df1=copy(x)
  df=df[,c('seqnames','start','end','pop')] %>% unique()
  df=df[,totsnps:=.N]
  
  df1=df1[ genomic_element=='Exons'][,c('seqnames','start','end','pop')] %>% unique()
  df1=df1[,totsnps_exons:=.N]
  
  df=df[df1,on=c('seqnames','start','end','pop'),nomatch=0]
  df=df[
    ,c('pop','totsnps','totsnps_exons')
  ] %>% unique()
return(df)
  
}

denisova_exons=exons(snps_annotated_proportion[[1]])
neandertal_exons=exons(snps_annotated_proportion[[2]])
png_exons=exons(snps_annotated_proportion[[3]])


test_exons=function(x){
  asnps=copy(x)
  nasnps=copy(png_exons)
  asnps=asnps[,
              asnp_nasnp_ratio:=totsnps/nasnps$totsnps][
    ,totnot_exons:=totsnps-totsnps_exons
  ]
  
  
  asnps=binom.test(x=c(asnps$totsnps_exons,nasnps$totsnps_exons),
                   p=asnps$asnp_nasnp_ratio,alternative = 'l')
  return(asnps)
}


test_exons(denisova_exons)
test_exons(neandertal_exons)

