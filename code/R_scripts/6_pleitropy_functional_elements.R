## look at pleitropic activity of snps across states
## expectations is that introgression is depleted from pleiotropic functional elements (Tellis et al 2019)
library(data.table);library(magrittr);library(dplyr)
library(GenomicRanges)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(R.utils);library(openxlsx)

options(width = 150)
numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

### load cells 
setwd('/data/projects/punim0586/dvespasiani/Files/')

table_dir='./Archaic_introgression_in_PNG/Results/Tables/'
plot_dir='./Archaic_introgression_in_PNG/Results/Plots/Chromatin_State/'


read_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            select = c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq','element_seqnames','element_start','element_end'))[
              ,freq_range:=ifelse(all_freq<0.05,'low','high')
              ]
    )
  pop_names=c('denisova','neandertal','png') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


snps_chrom_states=read_states('./Archaic_introgression_in_PNG/Chromatin_states/SNPs_chromHMM_annotated/new_set/')

chrom_state_levels=copy(snps_chrom_states[[1]])[,chrom_state] %>% unique() %>% stringr::str_sort( numeric = TRUE)


pleiotropy=function(x,freq){
  
  
  calculate_pleiotropy=function(x){
    x=x[,pleiotropy:=.N,by=.(seqnames,start,end,chrom_state)][,c('pleiotropy','chrom_state')] %>%
      group_by(`chrom_state`) %>% 
      mutate('median'=median(pleiotropy)) %>% 
      mutate('lowq'=quantile(pleiotropy,0.25)) %>% 
      mutate('upq'=quantile(pleiotropy,0.75))%>% as.data.table()
    
  }
  
  
  df=copy(x)
  
  if(freq=='none'){
  df1=copy(df)
  df1=lapply(df1,function(x)x=x[,c(1:4,6)] %>% unique())
  df1=lapply(df1,function(x)calculate_pleiotropy(x))
  df1=Map(mutate,df1,'pop'=names(df1)) %>% lapply(function(x)setDT(x))
                                                  
  return(df1)
  
  } else {
    df2=copy(df)
    df2=lapply(df2,function(x)x=x[freq_range%in%freq][,c(1:4,6)] %>% unique())
    df2=lapply(df2,function(x)calculate_pleiotropy(x))
    df2=Map(mutate,df2,'pop'=names(df2)) %>% lapply(function(x)setDT(x))
    
  return(df2)
  }

}
  
snps_all_freq_pleiotropy=pleiotropy(snps_chrom_states,'none')
snps_highfreq_pleiotropy=pleiotropy(snps_chrom_states,'high')

pvalues=function(x){
  df=copy(x) %>% rbindlist()
  df=df %>% split(as.factor(df$chrom_state))
  
  two_tailed_ttest=function(x){
    df=copy(x)
    deni=copy(df)[pop%in%'denisova']
    nean=copy(df)[pop%in%'neandertal']
    png=copy(df)[pop%in%'png']
    
    perform_test=function(x,y){
      x=x[,p:=t.test(x$pleiotropy,y$pleiotropy)$p.val
          ][
            ,statistic:=t.test(x$pleiotropy,y$pleiotropy)$stat
            ][
              ,df:=t.test(x$pleiotropy,y$pleiotropy)$parameter
              ][
                ,mean_aSNPs:=t.test(x$pleiotropy,y$pleiotropy)$estimate[[1]]
                ][
                  ,mean_naSNPs:=t.test(x$pleiotropy,y$pleiotropy)$estimate[[2]]
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
  
  df=lapply(df,function(x)two_tailed_ttest(x)) %>% rbindlist()
  df=df[,c('chrom_state','pop','mean_aSNPs','mean_naSNPs','statistic','df','p','p.adj','p.signif')]%>%unique()
  
  return(df)
}

pvalues_allfreq=pvalues(snps_all_freq_pleiotropy)
pvalues_highfreq=pvalues(snps_highfreq_pleiotropy)

pval=list(pvalues_allfreq,pvalues_highfreq)
pval=lapply(pval,function(x)x=x[order(factor(x$chrom_state,levels=chrom_state_levels))])

names(pval)=c('all frequencies','common-to-high frequency')
write.xlsx(pval,paste(table_dir,'Supp_Table_4_pleiotropy_pvalues.xlsx',sep=''))

## this is to plot values (bit messy with the first one but correct)
second_pleiotropy=function(x,freq){
  
  calculate_pleiotropy=function(x){
    df=copy(x)
    df=lapply(df,function(x)x=x[,pleiotropy:=.N,by=.(seqnames,start,end,chrom_state)][
      ,totsnps_chromstate:=.N,by=.(chrom_state)
      ][,totsnps_pleiotropy_chromstate:=.N,by=.(chrom_state,pleiotropy)][
        ,fraction:=totsnps_pleiotropy_chromstate/totsnps_chromstate
        ][,c('chrom_state','pleiotropy','fraction')] %>% unique())
    
    
    df=lapply(df,function(x)x=x[
      ,pleiotropy_ranges:= ifelse(pleiotropy==1 ,'1',
                                  ifelse(pleiotropy>=2 &pleiotropy<5,'2-4',
                                         ifelse(pleiotropy>=5 &pleiotropy<9,'5-8',
                                                ifelse(pleiotropy>=9 &pleiotropy<13,'9-12','13-18'))))][
                                                  ,c('pleiotropy_ranges','chrom_state','fraction')]%>%
        group_by(`chrom_state`,`pleiotropy_ranges`) %>% 
        mutate('median'=median(fraction)) %>% 
        mutate('lowq'=quantile(fraction,0.25)) %>% 
        mutate('upq'=quantile(fraction,0.75))%>% as.data.table())
    return(df)
  }
  
  df=copy(x)
  df=lapply(df,function(x)x=x[,c('element_seqnames','element_start','element_end','cell_type'):=NULL] %>% unique())
  if(freq=='none'){
    df1=copy(df)
    df1=calculate_pleiotropy(df1)
    df1=Map(mutate,df1,'pop'=names(df1)) %>% rbindlist() %>% unique()
    return(df1)
  }else{
    df2=copy(df)
    df2=lapply(df2,function(x)x=x[freq_range%in%freq])%>% calculate_pleiotropy()
    df2=Map(mutate,df2,'pop'=names(df2)) %>% rbindlist() %>% unique()
    return(df2)
  }
 
}

second_snps_all_freq_pleiotropy=second_pleiotropy(snps_chrom_states,'none')
second_snps_high_freq_pleiotropy=second_pleiotropy(snps_chrom_states,'high')


second_pleiotropy_plot=function(df){
    x=copy(df)
    x=x[,chrom_state:=gsub(".*_","",chrom_state)]

    x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
                                                  'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
    x$pleiotropy_ranges=factor(x$pleiotropy_ranges,levels = c('1','2-4','5-8','9-12','13-18'))
    my_palette_pop=c('#C99E10', # denisova
                     '#9B4F0F', #neandertal
                     '#1E656D' #png
                     )

    names(my_palette_pop)= levels(as.factor(x$pop))
    dotcol=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neanderthal','PNG'))
    barcol=scale_color_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neanderthal','PNG'))

    plot=ggplot(x,aes(x=pleiotropy_ranges,y=median,fill=pop))+
      geom_line(aes(group=pop,col=pop),size=0.3,position =position_dodge(width = 0.5))+
      geom_errorbar(aes(ymin=lowq, ymax=upq,col=pop),position=position_dodge(width = 0.5),width=0,size=1.2)+
          geom_point(size=3,position=position_dodge(width =0.5))+
          geom_point(shape = 21,colour = "black",position=position_dodge(width = 0.51),stroke = 0.25,size=3.97)+
      dotcol+barcol+
     xlab(' ')+ylab('\n Fraction SNPs for pleiotropy range \n ')+
      # geom_text(aes(x=pop, y=15.5, label=p.signif), col='black', size=3)+
     facet_wrap(chrom_state~., ncol = 3,scales = 'free_y')+
      theme(strip.text.x = element_text(),
            strip.text.y = element_text(hjust = 0.5),
            strip.background = element_rect(color = 'black', linetype = 'solid'),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
           panel.background =element_rect(fill = 'white', size = 0.25,colour = 'black'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.text = element_text(),
            legend.title = element_text(),
            legend.position = 'bottom',
           axis.text.y = element_text(),
            axis.title.y = element_text(hjust=0.5),
            axis.text=element_text(),
            axis.line = element_blank())
    return(plot)

  }

pdf(paste(plot_dir,'pleiotropy/pleiotropy_allfreq_plot.pdf',sep=''),width = 7,height = 7)
second_pleiotropy_plot(second_snps_all_freq_pleiotropy)
dev.off()

pdf(paste(plot_dir,'pleiotropy/pleiotropy_highfreq_plot.pdf',sep=''),width = 7,height = 7)
second_pleiotropy_plot(second_snps_high_freq_pleiotropy)
dev.off()
