## first figure of novel vs known variants and predicted impact ##
library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(R.utils)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)

setDTthreads(10)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

read_vep=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      y=fread(y,header = T,sep=' ')[,freq_range:=ifelse(all_frequency<0.05,'low','high')] 
    )
  x=x[c(2:4)]
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

ooa_snps=read_vep('./OoA_snps/')

## variant distribution across genomic elements
genomic_elements_perfreqbins=function(x){
    x=x[,genomic_element:=plyr::revalue(x$genomic_element,c("3_prime_UTR_variant"="3UTR",
                                                                     'stop_lost'='3UTR',
                                                                     "5_prime_UTR_variant" ="5UTR",
                                                                     'start_lost'='5UTR',
                                                                     "upstream_gene_variant"='Regulatory',
                                                                     "downstream_gene_variant"='Regulatory',
                                                                     "regulatory_region_variant"='Regulatory',
                                                                     "TF_binding_site_variant"='Regulatory',
                                                                     'mature_miRNA_variant'='Regulatory',
                                                                     "intergenic_variant"='Intergenic',
                                                                     "intron_variant"='Intron',
                                                                     "missense_variant" ='Exon',
                                                                     'synonymous_variant'='Exon',
                                                                     'stop_gained'='Exon',
                                                                     'stop_retained_variant'='Exon',
                                                                     'incomplete_terminal_codon_variant'='Exon',
                                                                     'coding_sequence_variant'='Exon',
                                                                     "non_coding_transcript_exon_variant"='Exon',
                                                                     "splice_acceptor_variant"='Splice',
                                                                     'splice_donor_variant'='Splice',
                                                                     "splice_region_variant"='Splice'))]
x=x[, c('seqnames','start','end','all_frequency','pop','genomic_element')][,all_frequency:=round(all_frequency,1)] %>% unique()

}

ooa_snps_genomicelement=copy(ooa_snps)
ooa_snps_genomicelement=lapply(ooa_snps_genomicelement,function(x)genomic_elements_perfreqbins(x))

fraction_totsnps_genomic_element=function(x){
  df=copy(x)[,all_frequency:=NULL] %>% unique()
  df=df[,totsnsp:=.N][,totsnps_genomicelement:=.N,by=.(genomic_element)][,fraction:=round((totsnps_genomicelement/totsnsp)*100,1)]
  df=df[,c('pop','fraction','genomic_element')] %>% unique()
  return(df)
}

snp_impact_per_element=lapply(ooa_snps_genomicelement,function(x)fraction_totsnps_genomic_element(x))


fraction_totsnps_genomic_element_per_frequency=function(x){
  df=copy(x)
  df=lapply(df,function(y)y=y[,'number_snps_perfreqbin':= .N, by=all_frequency] [
      ,'numbsnps_pergenomicelement_perfreqbin':= .N, by=.(all_frequency,genomic_element)][
        ,'fraction_pergenomicelement_perfreqbin':= numbsnps_pergenomicelement_perfreqbin/number_snps_perfreqbin
        ][
          ,c(4:6,9)]) %>% rbindlist()
  df=df[,pop:= plyr::revalue(df$pop,c('denisova_png'='Denisova',
                                                                        'neandertal_png'='Neandertal',
                                                                        'non_archaics_png'='Papuans'))][
                                                                          ,c('fraction_pergenomicelement_perfreqbin','all_frequency','genomic_element','pop')
                                                                          ] %>% unique()
  
  return(df)
}


snp_impact_per_element_perfreq=fraction_totsnps_genomic_element_per_frequency(ooa_snps_genomicelement)

my_palette_genomicelement=c('darkorchid1', # 3utr
             'cadetblue2', #5utr
             'tan3', #exons
             'thistle3', #intergenic
             'olivedrab3',
             'lightgoldenrod', # Regulatory
             'sienna4' #splice
             
)

names(my_palette_genomicelement)=levels(as.factor(x$genomic_element))
colScale_genomicelement=scale_fill_manual(name= " ",values = my_palette_genomicelement,
                                          labels = c("3UTR",
                                                     "5UTR",
                                                     "Exon",
                                                     "Intergenic",
                                                     "Intron",
                                                     "Regulatory",
                                                     "Splice"))

genomicelement_plot=ggplot(snp_impact_per_element_perfreq, aes(fill=genomic_element, y=fraction_pergenomicelement_perfreqbin, x=all_frequency)) + 
  geom_bar(stat="identity",col='black')+colScale_genomicelement+
  ylab('\n Fraction of SNPs \n')+ xlab('\n Allele frequency in PNG \n')+
  facet_wrap(pop~.,ncol=1)+
  theme(strip.text.y = element_text(),
        strip.text.x = element_text(),
       strip.background = element_blank(),
       strip.background.y = element_blank(),
       strip.background.x =element_blank(),
       panel.spacing=unit(1, "lines"),
      panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
     panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       axis.line = element_line(color = "black", size =0, linetype = "solid"))

pdf('/home/dvespasiani/snps_genomic_element_barplot.pdf',width = 7,height =7)
genomicelement_plot
dev.off()

