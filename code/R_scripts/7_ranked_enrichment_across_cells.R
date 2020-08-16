## rank enrichment per each state across cells and see where immune cells end up being
library(data.table);library(magrittr);library(dplyr)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(openxlsx)

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/Chromatin_states/')

# plot_output_dir='~/Archaic_introgression_in_PNG/ranked_enrichments/'
plot_dir='../Results/Plots/Chromatin_State/relative_enrichment/'
table_dir='../Results/Tables/'

assign_names=function(x){
  pop_names=c('denisova','neandertal') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

read_states=function(x){
  dir=paste0('./simplified_set/new_set/',x,sep='')
  
  y=as.character(list.files(dir,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T) %>% unique())
  y=lapply(y,function(z)
    z=z[
      ,cell_line:=plyr::revalue(z$cell_line,c("IMR90_fetal_lung_fibroblast"="IMR90",
                                              'BCells'='B cells',
                                              'TCells'='T cells',
                                              'Other_cells'= 'Other cells',
                                              'Smooth_muscle'='Smooth muscle',
                                              'ES_cells'='ES cells',
                                              'ES_derived_cells'='ES derived cells'))
      ])
  y=assign_names(y) 
  return(y)
}

allfreq_snps_chromstates=read_states('all_freq')
highfreq_snps_chromstates=read_states('high_freq')


chrom_state_levels=copy(allfreq_snps_chromstates[[1]])[,chrom_state] %>% unique() %>% stringr::str_sort( numeric = TRUE)

cell_line_enrichment=function(x){
  df=copy(x)
  df=lapply(df,function(a)a=a[
    ,numbsnps_perstate_percelltype:=NULL
    ][
      ,c('chrom_state','cell_type','cell_line','asnps_nasnps_ratio','mean_asnps_nasnps_ratio_chromstate')
      ][,enrichment:=log2(asnps_nasnps_ratio/mean_asnps_nasnps_ratio_chromstate)
        ][
          ,mean:=mean(enrichment),by=.(chrom_state,cell_line)][
            ,mean:=ifelse(is.na(mean),0,mean)
          ][
            ,sd:=sd(enrichment),by=.(chrom_state,cell_line)
            ][
              ,sd:=ifelse(is.na(sd),as.numeric(0),sd)] %>% unique())
  
  return(df)
  
}

asnp_allfreq_enrichment=cell_line_enrichment(allfreq_snps_chromstates)
asnp_highfreq_enrichment=cell_line_enrichment(highfreq_snps_chromstates)


## t-test on log2(enrichment)
t_test_cell_line=function(y){
  perform_test=function(x){
    df=copy(x)
    Tcells=copy(df)[cell_line%in%'T cells']
    Allcells=copy(df)[!cell_line%in%'T cells']
    df=df[
      ,p:=t.test(Tcells$enrichment,Allcells$enrichment,alternative = 'g')$p.value
      ][
        ,statistic:=t.test(Tcells$enrichment,Allcells$enrichment,alternative = 'g')$statistic	
        ][
          ,df:=t.test(Tcells$enrichment,Allcells$enrichment,alternative = 'g')$parameter
          ]
  }
  
  z=copy(y)
  z=z %>% split(as.factor(z$chrom_state)) %>% lapply(function(a)perform_test(a)) %>% rbindlist()
}

simplified=function(x){
  z=copy(x)
  z=lapply(z,function(a)a=a[,c(2:8):=NULL] %>% unique())
  z=Map(mutate,z,'pop'=names(z)) %>% rbindlist()
}


adjust_pvalues=function(x){
  df=copy(x)
  pvals_df=copy(df)
  pvals_df=pvals_df$p
  pvals_df_adjusted=p.adjust(pvals_df,'fdr')%>%as.data.table() %>% setnames('p.adj')
  pvals_df_adjusted=pvals_df_adjusted[
    ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                       ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                              ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                     ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
    ]
  df_final=cbind(df,pvals_df_adjusted)
  return(df_final)
}



asnp_allfreq_enrichment_pval=lapply(asnp_allfreq_enrichment,function(x)t_test_cell_line(x)%>% adjust_pvalues()) 
asnp_allfreq_enrichment_pval=simplified(asnp_allfreq_enrichment_pval)

asnp_highfreq_enrichment_pval=lapply(asnp_highfreq_enrichment,function(x)t_test_cell_line(x)%>% adjust_pvalues())
asnp_highfreq_enrichment_pval=simplified(asnp_highfreq_enrichment_pval)

pval=list(asnp_allfreq_enrichment_pval,asnp_highfreq_enrichment_pval)
pval=lapply(pval,function(x)x=x[order(factor(x$chrom_state,levels=chrom_state_levels))])
write.xlsx(pval,paste(table_dir,'Supp_Table_5_immune_enrichment.xlsx',sep=''))


## plot enrichment as dot plots
## colors
nihroadmap_colors=copy(asnp_allfreq_enrichment[[1]])[
  ,chrom_state:=gsub(".*_","",chrom_state)
  ][
    ,col:=plyr::revalue(`chrom_state`,c('TssA'='#FF0000','TssAFlnk'='#FF6E00','TxFlnk'='#32CD32','Tx'='#008000',
                                        'TxWk'='#006400','EnhG'='#C2E105','Enh'='#FFFF00',
                                        'ZNF/Rpts'='#66CDAA','Het'='#8A91D0','TssBiv'='#CD5C5C','BivFlnk'='#E9967A',
                                        'EnhBiv'='#BDB76B','ReprPC'='#3A3838','ReprPCWk'='#808080',
                                        'Quies'='#DCDCDC'))
    ][
      ,c('chrom_state','col')
      ] %>% unique()

my_palette_states=nihroadmap_colors$col
names(my_palette_states)=nihroadmap_colors$chrom_state

colScale=scale_fill_manual(name= " ",values = my_palette_states,labels = names(my_palette_states))
dotcol=scale_color_manual(name= " ",values = my_palette_states,labels = names(my_palette_states))


prepare_data=function(df){
  x=copy(df)
  x=x[,c(1,3,7,8)]%>% unique()
  x=x[
    ,chrom_state:=gsub(".*_","",chrom_state)
    ][
      ,chrom_state:=factor(chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
                                      'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
    ][
      ,cell_line:=factor(cell_line,levels = c('Adipose','B cells','Brain','Digestive','Epithelial','ES cells','ES derived cells','Heart',
                                                'IMR90','iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other cells','Smooth muscle',
                                                'T cells','Thymus'))
    ]
  # x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
  #                                               'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
  # 
  # x$cell_line=factor(x$cell_line,levels = c('Adipose','B cells','Brain','Digestive','Epithelial','ES cells','ES derived cells','Heart',
  #                                           'IMR90','iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other cells','Smooth muscle',
  #                                           'T cells','Thymus'))
  
}

asnps_allfreq_cellenrich=lapply(asnp_allfreq_enrichment,function(x)prepare_data(x)) 
asnps_highfreq_cellenrich=lapply(asnp_highfreq_enrichment,function(x)prepare_data(x))


enrichment_plot=function(df){
   plot1=ggplot(df,aes(x=cell_line,y=mean,col=chrom_state))+
      geom_line(aes(group=1),color='black',size=0.3)+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position='dodge',width=0,size=1.2)+
      geom_point(size=4.5)+
      geom_point(shape = 21,colour = "black",stroke = 0.25,size=4.97)+
      dotcol+
      colScale+
      xlab('\n  \n')+ylab('\n Log2 aSNPs enrichment \n')+
      facet_wrap(chrom_state~., ncol = 3,scales = 'free_y')+
      theme(strip.text.y = element_blank(),
            strip.background = element_blank(),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
            panel.spacing.x = unit(0.1, "lines"),
            panel.spacing=unit(0.01, "lines"),
            panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
            axis.text.x = element_text(angle=50,vjust = 0.9,hjust = 1),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_blank())
   return(plot1)
}

  
allstates_allfreq_enrichment_plot=lapply(asnps_allfreq_cellenrich,function(x)enrichment_plot(x))
allstates_highfreq_enrichment_plot=lapply(asnps_highfreq_cellenrich,function(x)enrichment_plot(x))


pdf(paste0(plot_dir,'deni_allfreq_allstates_enrichment.pdf',sep=''),width=10,height = 10)
allstates_allfreq_enrichment_plot[[1]]
dev.off()

pdf(paste0(plot_dir,'nean_allfreq_allstates_enrichment.pdf',sep=''),width=10,height = 10)
allstates_allfreq_enrichment_plot[[2]]
dev.off()

pdf(paste0(plot_dir,'deni_highfreq_allstates_enrichment.pdf',sep=''),width=10,height = 10)
allstates_highfreq_enrichment_plot[[1]]
dev.off()

pdf(paste0(plot_dir,'nean_highfreq_allstates_enrichment.pdf',sep=''),width=10,height = 10)
allstates_highfreq_enrichment_plot[[2]]
dev.off()


cool_states=copy(asnp_highfreq_enrichment_pval)[!p.signif %in% ' ' & pop=='denisova'][,chrom_state:=gsub(".*_","",chrom_state)]
cool_states=cool_states$chrom_state %>% unique()


deni_relevant_states=copy(asnps_highfreq_cellenrich[[1]])[chrom_state%in%cool_states]

pdf(paste0(plot_dir,'deni_highfreq_relevantstates_ranked_enrichment.pdf',sep=''),width=5,height = 15)
ggplot(deni_relevant_states,aes(x=cell_line,y=mean,col=chrom_state))+
  geom_line(aes(group=1),color='black',size=0.3)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position='dodge',width=0,size=1.2)+
  geom_point(size=4.5)+
  geom_point(shape = 21,colour = "black",stroke = 0.25,size=4.97)+
  dotcol+
  colScale+
  xlab('\n  \n')+ylab('\n Log2 aSNPs enrichment \n')+
  facet_wrap(chrom_state~., ncol = 1,scales = 'free_y')+
  theme(strip.text.y = element_blank(),
        strip.background = element_blank(),
        strip.background.y = element_blank(),
        strip.background.x =element_blank(),
        panel.spacing.x = unit(0.1, "lines"),
        panel.spacing=unit(0.01, "lines"),
        panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
        axis.text.x = element_text(angle=50,vjust = 0.9,hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_blank())

dev.off()

