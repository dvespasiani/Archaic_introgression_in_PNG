## motifbreak results
## % of TFBS SNPs
## genomic location
library(dplyr); library(data.table)
library(magrittr); 
library(purrr)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(annotatr)

## tfbs snp results ###
motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
plot_dir='./Results/Plots/Motifbreak/'

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

tfbs=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)fread(y,sep=' ',header = T))
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  
}

snps_tfbs=tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/new_set/combined')

# lapply(snps_tfbs,function(x)x=x[,c('seqnames','start','end')] %>% unique() %>% nrow())
## 39269 deni
## 22457 nean
## 294031 png
## annotate these snpa and look where they occurr
annotation_function = function(x){
annots = c('hg19_basicgenes', 'hg19_genes_intergenic')
annotations = build_annotations(genome = 'hg19', annotations = annots)

df=copy(x)
df=lapply(df,function(y)
  y=y[,c('seqnames','start','end')] %>% unique() %>%
    makeGRangesFromDataFrame())

df=lapply(df,function(y)  
  y=annotate_regions(regions =y,
                     annotations = annotations,
                     ignore.strand = T, quiet = F) %>% as.data.table())
df=lapply(df,function(y)y=y[,c('seqnames','start','end','annot.type')] %>% unique())
return(df)
}

tfbs_snps_annotated=annotation_function(snps_tfbs)


genome_location=function(x){
  df=copy(x)
  df=lapply(df,function(x)
    x=x[,tot_snps:=.N
        ][
          ,tot_snps_per_genome:=.N,by=.(annot.type)
          ][
            ,fraction:=tot_snps_per_genome/tot_snps
            ][
              ,c('fraction','annot.type')
              ] %>% unique() )
  df=Map(mutate,df,'pop'=names(df)) %>% rbindlist()
  
  df=df[
    ,genomic_element:=plyr::revalue(annot.type,c('hg19_genes_exons'='Exons',
                                                 'hg19_genes_introns'='Introns',
                                                 'hg19_genes_intergenic'='Intergenic',
                                                 'hg19_genes_promoters'='Promoters',
                                                 'hg19_genes_1to5kb'='Enhancers',
                                                 'hg19_genes_3UTRs'='3UTRs',
                                                 'hg19_genes_5UTRs'='5UTRs'))
    ][
      ,annot.type:=NULL
      ]
  return(df)
  
}
tfbs_genome_location=genome_location(tfbs_snps_annotated)

mean(tfbs_genome_location[genomic_element=='Enhancers']$fraction)*100 ## 6.00 %
mean(tfbs_genome_location[genomic_element=='Promoters']$fraction)*100 ## 1.64 %

#barplot tfbs location
tfbs_location_plot=ggplot(tfbs_genome_location,
                          aes(x=reorder(genomic_element,-fraction),
                              y=fraction,fill=pop,col=genomic_element))+
  geom_bar(stat="identity",position ='dodge')+
  scale_fill_manual(values=c('#C99E10','#9B4F0F','#1E656D'),name='Population',
                    labels=c('Denisova','Neandertal','PNG'))+ 
  scale_color_manual(values=c(rep('black',7)))+
  ylab('\n Proportion TFBS SNPs \n')+
  xlab('\n  \n')+
  theme(panel.background =element_rect(fill = 'white', colour = 'white'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(),
        legend.title = element_text(),
        legend.margin = margin(c(0.5, 2, 8, 25)),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text.x = element_text(),
        axis.text.y = element_text(),
        axis.title.y = element_text(hjust=0.5),
        axis.text=element_text(),
        axis.line = element_line(color = "black",size = 0.5, linetype = "solid"))

pdf(paste(plot_dir,'tfbs_snps_location_plot.pdf',sep=''),width = 8,height = 5)
tfbs_location_plot
dev.off()