## look at pleitropic activity of snps across states
## expectations is that introgression is depleted from pleiotropic functional elements (Tellis et al 2019)]
library(data.table);library(magrittr);library(dplyr)
library(GenomicRanges)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(R.utils);library(openxlsx)

setDTthreads(8)

### load cells 
setwd('/data/projects/punim0586/dvespasiani/Files/')
keys=c('seqnames','start','end')

assign_names=function(x,list_df){
  x=as.character(list.files(x,recursive=F,full.names=F))
  names(list_df)=x
  for (i in seq_along(list_df)){
    assign(x[i],list_df[[i]],.GlobalEnv)}
  return(list_df)
}

read_cells=function(x){
  df=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y) 
      fread(y,sep = ' ',header=T) %>% 
        makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
        as.data.table()
    )
  df=assign_names(x,df)
  df=rbindlist(df)
  df=df[
    ,c('width','strand'):=NULL
    ][
    ,pleiotropy:=.N, by=.(seqnames,start,chrom_state)
    ]
  
  return(df)
}

cells=read_cells('./Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states')


read_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            select = c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq','element_seqnames','element_start','element_end'))[
              ,freq_range:=ifelse(all_freq<0.05,'low','high')
              ]
    )
  # x=x[c(2:4)]
  pop_names=c('denisova','neandertal','png') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


snps_chrom_states=read_states('./PNG/Chromatin_states/SNPs_chromHMM_annotated/new_set/')

snps_chrom_states=lapply(snps_chrom_states,function(x)
  x=inner_join(x,cells,by=c('seqnames','element_start'='start','element_end'='end','chrom_state','cell_type')) %>% as.data.table())

pleiotropy=function(y,split){
 pleiotropy_function=function(z){
   df=copy(z)
  df=lapply(df,function(x)
    x=x[
      ,c('element_seqnames','element_start','element_end','cell_type','cell_line'):=NULL
      ] %>% unique())
  df=lapply(df,function(x)
    x=x %>% group_by(`seqnames`,`start`,`end`,`chrom_state`) %>%
      summarise(pleiotropy=sum(pleiotropy)) %>%
      as.data.table())
  
  df=lapply(df,function(x)
    x=x[
      ,pleiotropy:=ifelse(pleiotropy>111,111,pleiotropy)
      ][
        ,pleiotropy_ranges:=ifelse(pleiotropy>=1 & pleiotropy<4,'1_3',
                                   ifelse(pleiotropy>3 &pleiotropy<11,'4_10',
                                          ifelse(pleiotropy>10 &pleiotropy<26,'11_25',
                                                 ifelse(pleiotropy>25 &pleiotropy<51,'26_50','51_111'))))
        ]
    )
  df=Map(mutate,df,'pop'=names(df)) %>% lapply(function(x)setDT(x))
 }
 
  if(split=='no'){
    y=pleiotropy_function(y)
    return(y)
  } else{ 
    y=lapply(y,function(x)x[freq_range%in%'high'])
    y=pleiotropy_function(y)
    return(y)
    }
}

pleiotropy_snps=pleiotropy(snps_chrom_states,'no')
pleiotropy_highfreq_snps=pleiotropy(snps_chrom_states,'high')


# test=copy(snps_chrom_states)
# test=lapply(test,function(x)
#   x=x[
#     ,c('element_seqnames','element_start','element_end','cell_type','cell_line'):=NULL
#   ] %>% unique())
# 
# test=lapply(test,function(x)
#   x=x %>% group_by(`seqnames`,`start`,`end`,`chrom_state`) %>%
#     summarise(pleiotropy=sum(pleiotropy)) %>%
#     as.data.table())
# 
# test=lapply(test,function(x)
#   x=x[
#     ,pleiotropy:=ifelse(pleiotropy>111,111,pleiotropy)
#     ][
#     ,pleiotropy_ranges:=ifelse(pleiotropy>=1 & pleiotropy<4,'1_3',
#                                ifelse(pleiotropy>3 &pleiotropy<11,'4_10',
#                                       ifelse(pleiotropy>10 &pleiotropy<26,'11_25',
#                                              ifelse(pleiotropy>25 &pleiotropy<51,'26_50','51_111'))))
#     ]
# )
# test=Map(mutate,test,'pop'=names(test)) %>% lapply(function(x)setDT(x))

pvalues=function(x){
  df=copy(x) %>% rbindlist()
  df=df %>% split(as.factor(df$chrom_state))
 
one_tailed_wilcoxon=function(x){
    df=copy(x)
    deni=copy(df)[pop%in%'denisova']
    nean=copy(df)[pop%in%'neandertal']
    png=copy(df)[pop%in%'png']
    
    perform_test=function(x,y){
      x=x[,p:=wilcox.test(x$pleiotropy,y$pleiotropy,exact = F,alternative = 'l')$p.val
          ][
            ,statistic:=wilcox.test(x$pleiotropy,y$pleiotropy,exact = F,alternative = 'l')$stat
            ][
              ,p.adj:=p.adjust(p,method = 'bonferroni')
              ][
                ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                                   ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                          ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                                 ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
                ]
    }
    
    deni=perform_test(deni,png)
    nean=perform_test(nean,png) 
    combined=rbind(deni,nean) 
    return(combined)
}

df=lapply(df,function(x)one_tailed_wilcoxon(x)) %>% rbindlist()
df=df[,c(1:3,5,6):=NULL]%>%unique()

  return(df)
}

pvalues_allfreq=pvalues(pleiotropy_snps)
pvalues_highfreq=pvalues(pleiotropy_highfreq_snps)

# pvalues=copy(test) %>% rbindlist()
# pvalues=pvalues %>% split(as.factor(pvalues$chrom_state))
# 
# one_tailed_wilcoxon=function(x){
#   df=copy(x)
#   deni=copy(df)[pop%in%'denisova']
#   nean=copy(df)[pop%in%'neandertal']
#   png=copy(df)[pop%in%'png']
#   
#   perform_test=function(x,y){
#     x=x[,p:=wilcox.test(x$pleiotropy,y$pleiotropy,exact = F,alternative = 'l')$p.val
#         ][
#           ,statistic:=wilcox.test(x$pleiotropy,y$pleiotropy,exact = F,alternative = 'l')$stat
#           ][
#             ,p.adj:=p.adjust(p,method = 'bonferroni')
#             ][
#               ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
#                                  ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
#                                         ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
#                                                ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
#               ]
#   }
#   
#   deni=perform_test(deni,png)
#   nean=perform_test(nean,png) 
#   combined=rbind(deni,nean) 
#   return(combined)
# }
# 

write.xlsx(pvalues_allfreq,'/home/dvespasiani/pvalue_tables/pleiotropy_pvalues_allfreq.xlsx')
write.xlsx(pvalues_highfreq,'/home/dvespasiani/pvalue_tables/pleiotropy_pvalues_highfreq.xlsx')


## plot
prepare_data=function(x,pvalues){
  df=copy(x)%>% rbindlist()
  df=df[,pleiotropy_ranges:=NULL] %>% unique() 
  df=df[,c('chrom_state','pleiotropy','pop')]
  
  pval=copy(pvalues)
  pval=pval[,c('pop','chrom_state','p.signif')] %>% unique()
  
  df2=full_join(df,pvalues,by=c('chrom_state','pop')) %>% setDT()
  df2=df2[
    ,chrom_state:=gsub('.*_','',chrom_state)
    ]
  df2[is.na(df2)]=' '
  return(df2)
}

allfreq=prepare_data(pleiotropy_snps,pvalues_allfreq)
highfreq=prepare_data(pleiotropy_highfreq_snps,pvalues_highfreq)


pleiotropy_violinplot=function(x){
  
  x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
                                                'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
  
my_palette_pop=c('#C99E10', # denisova
                 '#9B4F0F', #neandertal
                 '#1E656D' #png
)

  names(my_palette_pop)= levels(as.factor(x$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))

  ggplot(x,aes(x=pop,y=pleiotropy,fill=pop))+
    colScale+
    xlab(' ')+ylab('\n Element pleiotropy \n ')+
      geom_violin(trim=T,scale = "width")+
      geom_boxplot(width=.04, position =  position_dodge(width = 0.2),outlier.shape = NA)+
    geom_text(data=x, aes(x=pop, y=(111)+0.0000001, label=p.signif), col='black', size=7)+
    # stat_compare_means(method = "wilcox.test",label = "p.signif",p.adjust.method = "BH",
    #                     ref.group = "png",hide.ns =T,
    #                     label.y =max(x$pleiotropy)+0.1,
    #                     size=12)+
    # scale_y_continuous(labels=c("0","40",'80','111'))+
    facet_wrap(chrom_state~., ncol = 3)+
    theme(strip.text.x = element_text(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          # panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black",size = 0, linetype = "solid"))

}

pdf('/home/dvespasiani/pleiotropy/pleiotropy_allfreq_plot.pdf',width = 10,height = 20)
pleiotropy_violinplot(allfreq)
dev.off()

pdf('/home/dvespasiani/pleiotropy/pleiotropy_highfreq_plot.pdf',width = 10,height = 20)
pleiotropy_violinplot(highfreq)
dev.off()

# 
# ## second pleiotropy mean + sd
# pleiotropy_dotplot=function(x){
#   x=copy(x)
#   x=x[,mean:=mean(pleiotropy),by=.(chrom_state)
#       ][
#         ,sd:=sd(pleiotropy),by=.(chrom_state)
#         ][
#           ,c(1,3,7,8,9)
#           ] %>% unique()
# 
  # x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
  #                                               'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))

#   
#   my_palette_pop=c('#C99E10', # denisova
#                                     '#9B4F0F', #neandertal
#                                     '#1E656D' #png
#                    )
# 
# names(my_palette_pop)= levels(as.factor(x$pop))
# colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
# 
#   ggplot(x,aes(x=pop,y=mean,fill=pop))+
#     geom_point()+
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0)+
#     geom_text(data=x, aes(x=pop, y=(111)+0.0000001, label=p.signif), col='black', size=7)+
#     # geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),size=0.3,position=position_dodge(0.05))+
#     colScale+
#     stat_summary(fun=mean, geom="line", aes(group=1)) +
#     xlab(' ')+ylab('\n Element pleiotropy \n ')+
#    facet_wrap(chrom_state~., ncol = 3)+
#     theme(strip.text.x = element_text(),
#           strip.text.y = element_text(hjust = 0.5),
#           strip.background = element_rect(color = 'black', linetype = 'solid'),
#           strip.background.y = element_blank(),
#           strip.background.x =element_blank(),
#           # panel.spacing=unit(1, "lines"),
#           panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           legend.text = element_text(),
#           legend.title = element_text(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.text.y = element_text(),
#           axis.title.y = element_text(hjust=0.5),
#           axis.text=element_text(),
#           axis.line = element_line(color = "black",size = 0, linetype = "solid"))
#   
# }
# 
# pdf('/home/dvespasiani/pleiotropy/pleiotropy_allfreq_dotplot.pdf',width = 10,height = 20)
# pleiotropy_dotplot(allfreq)
# dev.off()
# 
# pdf('/home/dvespasiani/pleiotropy/pleiotropy_highfreq_dotplot.pdf',width = 10,height = 20)
# pleiotropy_dotplot(highfreq)
# dev.off()
