## first figure of novel vs known variants and predicted impact ##
library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(R.utils)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr);
library(ggExtra)

setDTthreads(10)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

read_vep=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      y=fread(y,header = T,sep=' ')[,freq_range:=ifelse(all_frequency<0.05,'low','high')] 
    )
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

ooa_snps=read_vep('./OoA_snps/')

## variant distribution across genomic elements
genomic_elements_perfreqbins=function(x){
    x=x[,genomic_element:=plyr::revalue(x$genomic_element,c("3_prime_UTR_variant"="UTR",
                                                                     'stop_lost'='UTR',
                                                                     "5_prime_UTR_variant" ="UTR",
                                                                     'start_lost'='UTR',
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
x=x[, c('seqnames','start','end','all_frequency','pop','genomic_element')] %>% unique()

}

ooa_snps_genomicelement=copy(ooa_snps)
ooa_snps_genomicelement=lapply(ooa_snps_genomicelement,function(x)genomic_elements_perfreqbins(x))


fraction_totsnps_genomic_element=function(x){
  df=copy(x)[,all_frequency:=NULL] %>% unique()
  df=df[,totsnsp:=.N][,totsnps_genomicelement:=.N,by=.(genomic_element)][,fraction:=round((totsnps_genomicelement/totsnsp)*100,1)]
  df=df[,c('pop','fraction','genomic_element')] %>% unique()
  return(df)
}

lapply(ooa_snps_genomicelement,function(x)fraction_totsnps_genomic_element(x))



propsnps_genomicelement_freqbin=function(x){
  df=copy(x)
  df=lapply(df,function(y)y=y[
    ,all_frequency:=round(all_frequency,1)
    ][
      ,'numb_snps_perfreqbin':= .N, by=all_frequency
      ][
        ,'numbsnps_genomicelement_freqbin':= .N, by=.(all_frequency,genomic_element)
        ][
          ,'prop_genomicelement_freqbin':= numbsnps_genomicelement_freqbin/numb_snps_perfreqbin
          ][
            ,log10numbsnps_perfreqbin:=log10(numb_snps_perfreqbin) ])%>% rbindlist()
    df=df[,pop:= plyr::revalue(df$pop,c('denisova'='Denisova',
                                      'neandertal'='Neandertal',
                                      'png'='Papuans'))][
                                        ,c('numb_snps_perfreqbin','log10numbsnps_perfreqbin',
                                           'prop_genomicelement_freqbin','all_frequency','genomic_element','pop')
                                        ] %>% unique()

  return(df)
}


snp_genomelemnt_perfreq=propsnps_genomicelement_freqbin(ooa_snps_genomicelement)
## add mock lines 
mock_deni=data.table(numb_snps_perfreqbin=c(0,0),
                     log10numbsnps_perfreqbin=c(0,0),
                     prop_genomicelement_freqbin=c(0,0),
                     all_frequency=c(0.9,1.0),
                     genomic_element=c('Intron','Intron'),
                     pop='Denisova')
mock_nean=data.table(numb_snps_perfreqbin=0,
                     log10numbsnps_perfreqbin=0,
                     prop_genomicelement_freqbin=0,
                     all_frequency=1.0,
                     genomic_element='Intron',
                     pop='Neandertal')
snp_genomelemnt_perfreq=rbind(snp_genomelemnt_perfreq,mock_nean,mock_deni)

vep_barplot=function(x,ancestry){
  
  df=copy(x)
  df=df[pop%in%ancestry]
  my_palette_genomicelement=c('tan3', #exons
                              'thistle3', #intergenic
                              'olivedrab3', #intron
                              'lightgoldenrod', # Regulatory
                              'sienna4', #splice
                              'darkorchid1' # 3utr
  )
  
  names(my_palette_genomicelement)=levels(as.factor(df$genomic_element))
  colScale_genomicelement=scale_fill_manual(name= " ",values = my_palette_genomicelement,
                                            labels = c("Exon",
                                                       "Intergenic",
                                                       "Intron",
                                                       "Regulatory",
                                                       "Splice",
                                                       'UTR'))
  
  df_plot=ggplot(df,aes(fill=genomic_element, y=prop_genomicelement_freqbin, x=all_frequency)) +
    geom_bar(stat="identity",col='black')+
    colScale_genomicelement+
    ylab('\n Fraction of SNPs \n')+ xlab('\n SNP frequency in PNG \n')+
    theme(legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"),
          panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black", size = 0, linetype = "solid"))
          
  return(df_plot)
  
}


pdf('~/Archaic_introgression_in_PNG/vep_plot/Denisova_vep_barplot.pdf',width = 7,height =7)
vep_barplot(snp_genomelemnt_perfreq,'Denisova')
dev.off()

pdf('~/Archaic_introgression_in_PNG/vep_plot/Neandertal_vep_barplot.pdf',width = 7,height =7)
vep_barplot(snp_genomelemnt_perfreq,'Neandertal')
dev.off()

pdf('~/Archaic_introgression_in_PNG/vep_plot/PNG_vep_barplot.pdf',width = 7,height =7)
vep_barplot(snp_genomelemnt_perfreq,'Papuans')
dev.off()


log10_numbsnps_plot=function(file,ancestry,pop_color){
  
  df=copy(file)
  df=df[pop%in%ancestry][,c('all_frequency','log10numbsnps_perfreqbin','pop')] %>% unique()
   
    df_plot=ggplot(df,aes(x=all_frequency,y=log10numbsnps_perfreqbin,fill=pop))+
      geom_bar(stat='identity',position=position_dodge())+
      scale_fill_manual(name= " ",values=pop_color)+
      xlab('')+
      ylab('\n Log10 number of aSNPs \n')+
      theme(
        legend.position = "none",
        legend.key = element_rect(fill = "white", colour = "black"),
        panel.background =element_rect(fill = 'white', colour = 'white',size = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"))
    
    
    return(df_plot)
  
}

pdf('~/Archaic_introgression_in_PNG/vep_plot/denisova_log10_numbsnps.pdf',width = 7,height =2)
log10_numbsnps_plot(snp_genomelemnt_perfreq,'Denisova','#C99E10')
dev.off()

pdf('~/Archaic_introgression_in_PNG/vep_plot/neandertal_log10_numbsnps.pdf',width = 7,height =2)
log10_numbsnps_plot(snp_genomelemnt_perfreq,'Neandertal','#9B4F0F')
dev.off()

pdf('~/Archaic_introgression_in_PNG/vep_plot/png_log10_numbsnps.pdf',width = 7,height =2)
log10_numbsnps_plot(snp_genomelemnt_perfreq,'Papuans','#1E656D')
dev.off()


## check if there is no statistical enrichment/depletion within exons

exons=function(x){
  df=copy(x)
  df=df[,totsnps:=.N][ genomic_element=='Exon'][,c(1:3,5,7)] %>% unique()
  df=df[
    ,numbsnps:=.N
  ]
  
}

denisova_exons=exons(ooa_snps_genomicelement[[1]])
neandertal_exons=exons(ooa_snps_genomicelement[[2]])
png_exons=exons(ooa_snps_genomicelement[[3]])


test_exons=function(archaic,nonarchaic,fraction_asnps_nasnp){
  x=copy(archaic)
  x=x[,c(1:5):=NULL][,fakecolumn:='fake'] %>% unique()
  y=copy(nonarchaic)
  y=y[,c(1:5):=NULL][,fakecolumn:='fake'] %>% unique()
  
  z=inner_join(x,y,by='fakecolumn') %>% as.data.table()
  z=z[
    ,totsnps_exons:=numbsnps.x+numbsnps.y
  ]
  
  z=binom.test(x=c(z$numbsnps.x,z$numbsnps.y),p=fraction_asnps_nasnp,alternative = 'l')
  return(z)
}


test_exons(denisova_exons,png_exons,143581/1112811)
test_exons(neandertal_exons,png_exons,85298/1112811)

# neandertal_exons=exons(ooa_snps_genomicelement[[2]])
# png_exons=exons(ooa_snps_genomicelement[[3]])




# ## log2 fold enrichment
# enrichment_per_genomic_element=function(x){
#   df=copy(x)
#   df=lapply(df,function(x)x=x[
#     ,all_frequency:=round(all_frequency,2)
#     ][
#       ,totsnsp_allfreq:=.N,by=.(all_frequency)
#       ][
#         ,totsnps_genomicelement_allfreq:=.N,by=.(genomic_element,all_frequency)
#         ][
#           ,fractionsnp_genomicelement_allfreq:=(totsnps_genomicelement_allfreq/totsnsp_allfreq)
#                                 ][
#                                   ,freq_bin:=round(all_frequency,1)
#                                 ])
# 
#   deni=copy(df[[1]])[,c('freq_bin','genomic_element','fractionsnp_genomicelement_allfreq','totsnsp_allfreq','pop')] %>% unique()
#   nean=copy(df[[2]])[,c('freq_bin','genomic_element','fractionsnp_genomicelement_allfreq','totsnsp_allfreq','pop')] %>% unique()
#   png=copy(df[[3]])[,c('freq_bin','genomic_element','fractionsnp_genomicelement_allfreq')] %>% unique()
# 
# 
#   arch_enrichment=function(x,y){
#     df=inner_join(x,y,by=c('freq_bin','genomic_element'))%>% as.data.table()
#     df=df[,log2_foldenrich_freqbin:=log2(fractionsnp_genomicelement_allfreq.x/fractionsnp_genomicelement_allfreq.y)
#         ][
#           ,c(3,6):=NULL
#           ] %>% unique()
#     df=df[
#       ,mean:=mean(log2_foldenrich_freqbin),by=.(freq_bin,genomic_element)
#           ][
#             ,sd:=sd(log2_foldenrich_freqbin),by=.(freq_bin,genomic_element)
#             ]
#     return(df)
#   }
# 
#   deni_enich=arch_enrichment(deni,png)
# 
# 
#  nean_enich=arch_enrichment(nean,png)
# 
# 
#  df_final=list(deni_enich,nean_enich)
#   return(df_final)
# }
# 
# aSNPs_enrichments=enrichment_per_genomic_element(ooa_snps_genomicelement)
# 
# 
# ## one tailed wilcox test on log2 fold enrichments (not adequate probably)
# ## do hypergeometric test
# stat_test=function(file){
#   df=copy(file)[
#   ,c('freq_bin','genomic_element','log2_foldenrich_freqbin')
# ] %>% unique()
# 
# df=df%>% split(as.factor(df$genomic_element)) %>%
#   lapply(function(x)x %>% split(as.factor(x$freq_bin)) %>%
#            lapply(function(y)
#              y=y[,p:=wilcox.test(y$log2_foldenrich_freqbin,alternative = 'g',mu = 0.0)$p.va][
#       ,statistic:=wilcox.test(y$log2_foldenrich_freqbin,alternative = 'g',mu = 0.0)$stat
#       ][
#         ,p.adj:=p.adjust(p,method = 'bonferroni')
#         ][
#           ,p.signif:= ifelse(`p.adj`<=0.05,'*',' ')
#                              # ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
#                              #        ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
#                              #               ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
#           ] ) %>% rbindlist()) %>% rbindlist()
# return(df)
# }
# 
# aSNPs_pvals=lapply(aSNPs_enrichments,function(x)stat_test(x))
# 
# aSNPs_enrichments_pvals=purrr::map2(aSNPs_enrichments,aSNPs_pvals,inner_join,by=c('freq_bin','genomic_element','log2_foldenrich_freqbin'))
# aSNPs_enrichments_pvals=lapply(aSNPs_enrichments_pvals,function(x)setDT(x))
# 
# 
# 
# foldenrich_plot=function(file,ancestry){
# 
#   plot_enrich=function(x){
#     my_palette_genomicelement=c('tan3', #exons
#                                 'thistle3', #intergenic
#                                 'olivedrab3', #intron
#                                 'lightgoldenrod', # Regulatory
#                                 'sienna4', #splice
#                                 'darkorchid1' # 3utr
#     )
# 
#     names(my_palette_genomicelement)=levels(as.factor(x$genomic_element))
#     colScale_genomicelement=scale_color_manual(name= " ",values = my_palette_genomicelement,
#                                                labels = c("Exon",
#                                                           "Intergenic",
#                                                           "Intron",
#                                                           "Regulatory",
#                                                           "Splice",
#                                                           'UTR'))
# 
#     x=x[,c('freq_bin','mean','sd','genomic_element','p.signif')] %>% unique()
# 
#     df_plot=ggplot(x,aes(x=freq_bin,y=mean,col=genomic_element))+
#       geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),position =position_dodge(width=0.07))+
#       geom_text(data=x,aes(x=freq_bin, y=mean+0.05, label=p.signif), col='black', size=5)+
#       geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.4)+
#       colScale_genomicelement+
#       facet_wrap(genomic_element~.,ncol = 1)+
#       xlab('\n aSNPs frequency bin \n')+
#       ylab('\n Log2 fold enrichment \n')+
#       theme(strip.text.y = element_blank(),
#             strip.background = element_blank(),
#             strip.background.y = element_blank(),
#             strip.background.x =element_blank(),
#             panel.spacing=unit(0.01, "lines"),
            # legend.position = "bottom",
            # legend.key = element_rect(fill = "white", colour = "black"),
            # panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
            # panel.grid.minor = element_blank(),
            # panel.grid.major = element_blank(),
            # axis.line = element_line(color = "black", size = 0, linetype = "solid"))
#     
#     return(df_plot)
#   }
# 
#   if(ancestry%in%'denisova'){
#     deni=copy(file[[1]])
#     deni_plot=plot_enrich(deni)
#     return(deni_plot)
#   }else{
#     nean=copy(file[[2]])
#     nean_plot=plot_enrich(nean)
#     return(nean_plot)
#   }
# }
# 
# 
# pdf('/home/dvespasiani/vep_plot/denisova_enrich.pdf',width = 5,height=10)
# foldenrich_plot(aSNPs_enrichments_pvals,'denisova')
# dev.off()
# 
# pdf('/home/dvespasiani/vep_plot/neandertal_enrich.pdf',width = 5,height=10)
# foldenrich_plot(aSNPs_enrichments_pvals,'neandertal')
# dev.off()
# 
# # '#C99E10', # denisova
# # '#9B4F0F', #neandertal
# log10_numbsnps_plot=function(file,ancestry){
# 
#   plot_log10=function(x,pop_color){
# 
#     df_snps=copy(x)[,c('freq_bin','totsnsp_allfreq','pop')] %>% unccique()
#     df_snps=df_snps[,totsnsp_allfreq:=log10(totsnsp_allfreq)]
# 
#     df_plot=ggplot(df_snps,aes(x=freq_bin,y=totsnsp_allfreq,fill=pop))+
#       geom_bar(stat='identity',position=position_dodge())+
#       scale_fill_manual(name= " ",values=pop_color)+
#       xlab('')+
#       ylab('\n Log10 number of aSNPs \n')+
#       theme(
#         legend.position = "none",
#         legend.key = element_rect(fill = "white", colour = "black"),
#         panel.background =element_rect(fill = 'white', colour = 'white',size = 0),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.line = element_line(color = "black", size = 0.5, linetype = "solid"))
# 
# 
#     return(df_plot)
#   }
#   if(ancestry%in%'denisova'){
#     deni=copy(file[[1]])
#     deni_plot=plot_log10(deni,'#C99E10')
#     return(deni_plot)
#   }else{
#     nean=copy(file[[2]])
#     nean_plot=plot_log10(nean,'#9B4F0F')
#     return(nean_plot)
#   }
# 
# }
# 
# 
# pdf('/home/dvespasiani/vep_plot/denisova_log10_numbsnps.pdf',width = 10,height=2)
# log10_numbsnps_plot(aSNPs_enrichments,'denisova')
# dev.off()
# 
# pdf('/home/dvespasiani/vep_plot/neandertal_log10_numbsnps.pdf',width = 10,height=2)
# log10_numbsnps_plot(aSNPs_enrichments,'neandertal')
# dev.off()



