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
combined_tfbs_outputdir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/combined/'
motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'


setDTthreads(8)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')


gtex_expr=fread('./GTEx/GTEx_Analysis_v8_gene_median_tmp/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',sep='\t',header = T)
gene_ids=gtex_expr[,c(1:2)]


read_tfbs_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
  fread(y,sep=' ',header = T)[
    ,c('start', 'end') :=NULL][
      ,start := snpPos][
        ,end := start +1][
          ,snpPos := NULL]
    )
  
  x=lapply(x,function(y)y=semi_join(y,gene_ids,by=c('geneSymbol'='Description')) %>% ## remove TFs absent from gtex (only D/OBOXs)
             as.data.table()) 
  
  pop_names=c('ambiguous','denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


hocomoco_tfbs=read_tfbs_snps(paste(motif_input_dir,'hocomoco',sep=''))
jaspar_tfbs=read_tfbs_snps(paste(motif_input_dir,'jaspar',sep=''))
# encode_tfbs=read_tfbs_snps(paste(motif_input_dir,'encode',sep=''))
# known_encode_tfbs=lapply(encode_tfbs,function(x)x=copy(x)[providerName%like%'_known_'])

combined=purrr::map2(hocomoco_tfbs,jaspar_tfbs,rbind)
# combined=purrr::map2(combined,hocomoco_tfbs,rbind)
combined=lapply(combined,function(x)x=unique(x))
# 
# # numb snps
# x=lapply(encode_tfbs,function(x)x=x[providerName%like%'REST'][,c('seqnames','start','end','providerName')] %>% unique())
# y=lapply(hocomoco_tfbs,function(x)x=x[providerName%like%'REST'][,c('seqnames','start','end','providerName')] %>% unique())
# 
# z=purrr::map2(x,y,semi_join,c('seqnames','start','end')) 
# z=lapply(z,function(a)a %>% unique()%>% nrow())
# 
# x=lapply(x,function(a)a=a[,c('seqnames','start','end')] %>% unique()%>% nrow())

# filenames=paste0(combined_tfbs_outputdir,names(combined),sep='')
# mapply(write.table, combined, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)


numb_tfbs_snps=function(x){lapply(x,function(y)y[,c('seqnames','start','end')]%>% unique() %>% nrow())}
## HOCOMOCO numb snps: 
## 2734 ambig 
## 18947 deni
## 10873 nean
## 58016 png

## JASPAR numb snps:
## 1753 ambig
## 12518 deni
## 7158 nean
## 37022 png 

## ENCODE numb snps:
## 5341 ambig
## 37961 deni
## 21531 nean
## 292689 png

## KNOWN ENCODE numb snps:
## 3745 ambig
## 26972 deni
## 15220 nean
## 206096 png

## TOTAL COMBINED 
## 5392 ambig
## 38730 deni
## 22080 nean
## 240389 png

## fraction of tfbs SNPs over total input snps
input_snps=as.character(list.files('./Motifbreak/Tx_and_CREs/snps_motifbreakr_format/non_splitted',recursive = F,full.names = T)) %>% 
  lapply(function(y)
    y=fread(y,sep='\t',header = F,drop='V4',col.names = c('seqnames','start','end'))
  )
## input number snps:
## 16034 Ambiguous
## 113674 Denisova
## 64794 Neandertal
## 878481 PNG

combined=lapply(combined,function(x)x=x[,c('width','motifPos','Refpvalue','Altpvalue','dataSource','seqMatch'):=NULL])

## annotate these snpa and look where they occurr
annotation_function = function(x){
annots = c('hg19_basicgenes', 'hg19_genes_intergenic')
annotations = build_annotations(genome = 'hg19', annotations = annots)

df=copy(x)
df=lapply(df,function(y)
  y=makeGRangesFromDataFrame(y,keep.extra.columns = T))

df=lapply(df,function(y)  
  y=annotate_regions(regions =y,
                     annotations = annotations,
                     ignore.strand = T, quiet = F) %>% as.data.table())
df=lapply(df,function(y)
          y=y[,c(1:3,6:8,17:20,26,27)] %>% unique())
return(df)
}

tfbs_snps_annotated=annotation_function(combined)


genome_location=function(x){
  df=copy(x)
  df=lapply(df,function(x)x=x[,c(1:3,12)] %>% unique())
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
tfbs_genome_location=tfbs_genome_location[!pop%in%'ambiguous']
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

pdf('/home/dvespasiani/tfbs_plots/tfbs_snps_location_plot.pdf',width = 8,height = 5)
tfbs_location_plot
dev.off()
