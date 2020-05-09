## rank enrichment per each state across cells and see where immune cells end up being
library(data.table);library(magrittr);library(dplyr)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(openxlsx)

setDTthreads(8)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/Chromatin_states/')

plot_output_dir='/home/dvespasiani/ranked_enrichments/'

assign_names=function(x){
  pop_names=c('denisova','neandertal','png') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

read_states=function(x){
  y=as.character(list.files(x,recursive = F,full.names = T)) %>% 
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

snps_chrom_states=read_states('./simplified_set/new_set')


enrichment=function(x){
  y=copy(x[[3]])[,c(4,5):=NULL]
  z=copy(x[c(1:2)]) %>% lapply(function(z)z=z[,c(4,5):=NULL]) 
  z=lapply(z,function(a)inner_join(a,y,by=c('chrom_state','cell_type','cell_line')) %>%as.data.table()) 
  z=lapply(z,function(a)a=a[
    ,enrichment :=fraction_snpdens_perstate_perepigenome.x/fraction_snpdens_perstate_perepigenome.y
    ][
      ,c('chrom_state','cell_type','cell_line','enrichment')
      ] %>% unique()
  )
  # z=assign_names(z)
  z=lapply(z,function(a)
    a=a[
      ,max_enrichment:=max(enrichment),by=.(chrom_state)
      ][
        ,ranked_enrichment:=enrichment/max_enrichment
        ][,enrichment:=log2(enrichment)]
  )
  return(z)
  
}


enrichment_highfreq=function(x){
  y=copy(x[[3]])[freq_range%in%'high'][,c(4,5):=NULL]
  z=copy(x[c(1:2)])
  z=lapply(z,function(a)a=a[freq_range%in%'high'][,c(4,5):=NULL]) 
  z=lapply(z,function(a)inner_join(a,y,by=c('chrom_state','cell_type','cell_line')) %>%as.data.table()) %>%
    lapply(function(a)a=a[
      ,enrichment :=fraction_snpdens_perstate_perepigenome_perfreqrange.x/fraction_snpdens_perstate_perepigenome_perfreqrange.y
      ][
        ,c('chrom_state','cell_type','cell_line','enrichment')
        ] %>% unique()
    )
  # z=assign_names(z)
  z=lapply(z,function(a)
    a=a[
      ,max_enrichment:=max(enrichment),by=.(chrom_state)
      ][
        ,ranked_enrichment:=enrichment/max_enrichment
        ][,enrichment:=log2(enrichment)]
  )
  return(z)
  
}

asnp_allfreq_enrichment=enrichment(snps_chrom_states)
asnp_highfreq_enrichment=enrichment_highfreq(snps_chrom_states)


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
          ][
            ,p.adj:=p.adjust(p,method = 'bonferroni')
            ][
              ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                                 ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                        ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                               ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
              ]
  }
  
  z=copy(y)
  z=z %>% split(as.factor(z$chrom_state)) %>% lapply(function(a)perform_test(a)) %>% rbindlist()
}

simplified=function(x){
  z=copy(x)
  z=lapply(z,function(a)a=a[,c(2:6):=NULL] %>% unique())
  z=Map(mutate,z,'pop'=names(z)) %>% rbindlist()
}

asnp_allfreq_enrichment_pval=lapply(asnp_allfreq_enrichment,function(x)t_test_cell_line(x))
asnp_allfreq_enrichment_pval=simplified(asnp_allfreq_enrichment_pval)

asnp_highfreq_enrichment_pval=lapply(asnp_highfreq_enrichment,function(x)t_test_cell_line(x))
asnp_highfreq_enrichment_pval=simplified(asnp_highfreq_enrichment_pval)

pval=list(asnp_allfreq_enrichment_pval,asnp_highfreq_enrichment_pval)

write.xlsx(pval,'/home/dvespasiani/pvalue_tables/Supp_Table_ranked_enrichment.xlsx')


### mean + sd 
# test=copy(asnp_allfreq_enrichment)
# 
# test=lapply(test,function(x)x=x[,mean:=mean(ranked_enrichment),by=.(chrom_state,cell_line)
#                                 ][,sd:=sd(ranked_enrichment),by=.(chrom_state,cell_line)][
#                                   ,c(1,3,7,8)
#                                 ][,sd:=ifelse(is.na(sd),0,sd)] %>% unique())

## plot enrichment as dot plots
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

ranked_enrichment_plot=function(x){
  x=copy(x)
  x=x[,mean:=mean(ranked_enrichment),by=.(chrom_state,cell_line)
      ][
        ,sd:=sd(ranked_enrichment),by=.(chrom_state,cell_line)
        ][
          ,c(1,3,7,8)
          ][
            ,sd:=ifelse(is.na(sd),0,sd)
            ][
              ,mean:=ifelse(is.na(mean),0,mean)] %>% unique()
  x=x[
    ,chrom_state:=gsub(".*_","",chrom_state)
    ]
  x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
                                                'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
  
  x$cell_line=factor(x$cell_line,levels = c('Adipose','B cells','Brain','Digestive','Epithelial','ES cells','ES derived cells','Heart',
                                            'IMR90','iPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other cells','Smooth muscle',
                                            'T cells','Thymus'))
  
  
  my_palette_states=nihroadmap_colors$col
  names(my_palette_states)=nihroadmap_colors$chrom_state
  
  colScale=scale_fill_manual(name= " ",values = my_palette_states,labels = names(my_palette_states))
  linecol=scale_color_manual(name= " ",values = my_palette_states,labels = names(my_palette_states))
  # linecol='black'
  
  ggplot(x,aes(x=cell_line,y=mean,col=chrom_state))+
    # geom_point()+
    # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(width=0.6),width=0)+
    geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),size=0.3)+
    linecol+
    colScale+
    geom_line(aes(group=1),color='black',size=0.3)+
    # stat_summary(fun=mean, geom="line", aes(group=1)) +
    xlab('\n  \n')+ylab('\n Relative enrichment \n')+
    ylim(0,1.01)+
    facet_wrap(chrom_state~., ncol = 1)+
    theme(strip.text.y = element_blank(),
          strip.background = element_blank(),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          panel.spacing=unit(0.01, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size = 0.5),
          axis.text.x = element_text(angle=50,vjust = 0.9,hjust = 1),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black", size =0, linetype = "solid"))
  
}
allstates_allfreq_enrichment_plot=lapply(asnp_allfreq_enrichment,function(x)ranked_enrichment_plot(x))
allstates_highfreq_enrichment_plot=lapply(asnp_highfreq_enrichment,function(x)ranked_enrichment_plot(x))

save_plots=function(plots,directory,plotname,width,height){
  lapply(names(plots),function(x)
    ggsave(filename=paste(directory,paste0(x,plotname,sep="")),limitsize = F,height = height,width = width,units = 'cm',
           plot=plots[[x]]))
  
}
save_plots(allstates_allfreq_enrichment_plot,plot_output_dir,'_allfreq_allstates_ranked_enrichment.pdf',10,50)
save_plots(allstates_highfreq_enrichment_plot,plot_output_dir,'_highfreq_allstates_ranked_enrichment.pdf',10,50)

## select cool states for the main
cool_states=c('7_Enh','6_EnhG','5_TxWk','4_Tx','8_ZNF/Rpts','2_TssAFlnk','14_ReprPCWk')

allfreq_asnp_enrichment_cool=copy(asnp_allfreq_enrichment) %>% lapply(function(x)x=x[chrom_state%in%cool_states])
# highfreq_asnp_enrichment_cool=copy(asnp_highfreq_enrichment) %>% lapply(function(x)x=x[chrom_state%in%cool_states])

allfreq_cool_states_enrichment_plot=lapply(allfreq_asnp_enrichment_cool,function(x)ranked_enrichment_plot(x))
# highfreqfreq_cool_states_enrichment_plot=lapply(highfreq_asnp_enrichment_cool,function(x)ranked_enrichment_plot(x))

save_plots(allfreq_cool_states_enrichment_plot,plot_output_dir,'_allfreq_relevant_states_ranked_enrichment.pdf',10,20)
# save_plots(highfreqfreq_cool_states_enrichment_plot,plot_output_dir,'_highfreq_activestates_ranked_enrichment.pdf',10,20)




