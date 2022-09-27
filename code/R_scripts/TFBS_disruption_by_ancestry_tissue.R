## check disruptive potential TFBS SNPs by ancestry
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(annotatr);library(data.table);library(dplyr)
library(openxlsx)
## tfbs snp results ###
motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

plot_dir='./Results/Plots/Motifbreak/'
table_dir='./Results/Tables/'
merging_keys=c('seqnames','start','end')

immune_cells='TCells'

cell_line_levels=c('Adipose','BCells','Brain','Digestive','Epithelial','ES_cells','ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
                   'IPSC','Mesenchymal','Muscle','Myosatellite','Neurospheres','Other_cells','Smooth_muscle','TCells','Thymus')

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')


read_tfbs_snps=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)fread(y,sep=' ',header = T))
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  
}

tfbs=read_tfbs_snps(paste(motif_input_dir,'new_set/combined/',sep=''))
lapply(tfbs,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

read_files=function(x){
  x=as.character(list.files(x,recursive =F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep = ' ',header = T)[
        ,pop:=NULL
        ] %>% 
        setnames(old=c('CHR','FROM','TO'),new = c('seqnames','start','end')))
  
  pop_names=c('denisova','neanderthal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  
}

original_snps=read_files('./Grouped_filtered_snps/new_set/') 

snps_tfbs=copy(tfbs)
snps_tfbs=purrr::map2(snps_tfbs,original_snps,inner_join,by=c('seqnames','start','end','REF','ALT')) %>% 
  lapply(function(x)setDT(x))

snps_tfbs=lapply(snps_tfbs,function(x)x=x[
  ,alleleDiff_adj:=ifelse((ancestry=='archaic' & POP_ARCH_ALT>POP_ARCH_REF),alleleDiff,
                          ifelse((ancestry=='non_archaic'& ANC=='0'),alleleDiff,-alleleDiff))
])
lapply(snps_tfbs,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

### look at effects tissue-wise
read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq'))[
        chrom_state%in%c('2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
        ][
          ,chrom_state:=NULL
          ][
            ,cell_line:=ifelse(cell_line=='iPSC','IPSC',cell_line)
          ]
    )
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') 

## combine tfbs with cell-type/state info
tfbs_cres=copy(snps_tfbs)
tfbs_cres=lapply(tfbs_cres,function(x)x=x[, .SD[which.max(abs(alleleDiff_adj))], by=.(seqnames,start,end,REF,ALT)])

lapply(snps_tfbs,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
lapply(tfbs_cres,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

tfbs_cres=purrr::map2(tfbs_cres,states,inner_join,by=merging_keys)
tfbs_cres=Map(mutate,tfbs_cres,'pop'=names(tfbs_cres))%>% 
  lapply(function(x)x=as.data.table(x)[
    ,c('seqnames','start','end','REF','ALT','cell_type','cell_line','all_freq','effect','alleleDiff_adj','pop')
    ] %>% unique())

lapply(tfbs_cres,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())


strong_highfreq_effect=copy(tfbs_cres)
strong_highfreq_effect=lapply(strong_highfreq_effect,function(x)x=x[,all_freq:=round(all_freq,1)][all_freq>=0.3][effect=='strong'])
lapply(strong_highfreq_effect,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

strong_highfreq_effect=rbindlist(strong_highfreq_effect)


### table p-values
list_tfbs=copy(strong_highfreq_effect)

stat_results=compare_means(
  alleleDiff_adj~pop,
  list_tfbs,
  method = "wilcox.test",
  paired = FALSE,
  group.by = 'cell_line',
  ref.group = 'png',
  p.adjust.method = "fdr"
) %>% as.data.table()

stat_results=stat_results[
  ,p.adj:=p.adjust(p,method = 'fdr'),by=.(group2)
  ][
    ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                       ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                              ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                     ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
    ][
      ,c('.y.','group1','p.format','method'):=NULL
      ] %>% setnames(old='group2',new='pop') %>% setorderv('cell_line',1)


get_median_effects=function(x,population){
  df=copy(x)[pop%in%population]
  df=df%>%group_by(cell_line) %>%
    summarize(int = median(alleleDiff_adj)) %>% 
    as.data.table() %>% setnames(old='int','median_allele_effect')
  df=df[,pop:=population][,median_allele_effect:=round(median_allele_effect,2)]
  df1=copy(stat_results)
  
  if(population=='png'){
    df=df[,pop:=NULL] %>% setnames(old='median_allele_effect',new='png_median_allele_effect')
    return(df)
  }else{
    df1=df1[df,on=c('cell_line','pop'),nomatch=0]
  }
  return(df1)
}

deni_disruption=get_median_effects(strong_highfreq_effect,'denisova')
nean_disruption=get_median_effects(strong_highfreq_effect,'neandertal')
png_disruption=get_median_effects(strong_highfreq_effect,'png')

stat_results_final=rbind(deni_disruption,nean_disruption)
stat_results_final=stat_results_final[png_disruption,on='cell_line',nomatch=0]

table_dir='./Results/Tables/'
write.xlsx(stat_results_final,paste(table_dir,'Supp_Table_9_TFs_effect_by_tissue.xlsx',sep=''))


## plot effect
plot_allele_effect=function(x,cells){

  df=copy(x)
  df$cell_line=factor(df$cell_line,levels=cell_line_levels)
  df=df[cell_line%in%cells]
  
  png=copy(df)[pop%in%'png']
  png=png%>%group_by(cell_line) %>%summarize(int = median(alleleDiff_adj))
  
  
  my_palette_pop=c(
    '#C99E10', # denisova
    '#9B4F0F', #neandertal
    '#1E656D'#png
  )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
  
  ggplot(df,aes(x=pop,y=alleleDiff_adj,fill=pop))+
    colScale+
    xlab(' ')+ylab('\n Delta PWM \n ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_text(data=stat_results_final, aes(x=pop, y=max(df$alleleDiff_adj)+0.1, label=p.signif), col='black', size=7)+
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey",size=0.5)+
    facet_wrap(cell_line~., ncol = 4)+
    geom_hline(data = png, aes(yintercept=int), linetype="dashed", color = "lightgrey",size=0.5)+
    theme(strip.text.x = element_text(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          legend.margin = margin(c(0.5, 2, 8, 25)),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black", size = 0, linetype = "solid"))

}

pdf(paste(plot_dir,'TFBS_effect_by_tissue.pdf',sep=''),width = 20,height = 20)
plot_allele_effect(strong_highfreq_effect,strong_highfreq_effect$cell_line)
dev.off()

## plot effect t cells (too lazy to change the above function to make it universal)
plot_allele_effect_tcells=function(x){
  
  df=copy(x)[cell_line%in%'TCells']
  png=copy(df)[pop%in%'png']
  png=png%>%group_by(cell_line) %>%summarize(int = median(alleleDiff_adj))
  
  stat_results_final_immune=stat_results_final[cell_line=='TCells']
  my_palette_pop=c(
    '#C99E10', # denisova
    '#9B4F0F', #neandertal
    '#1E656D'#png
  )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
  
  ggplot(df,aes(x=pop,y=alleleDiff_adj,fill=pop))+
    colScale+
    xlab(' ')+ylab('\n Delta PWM \n ')+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2)+
    geom_text(data=stat_results_final_immune, aes(x=pop, y=max(df$alleleDiff_adj)+0.1, label=p.signif), col='black', size=7)+
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey",size=0.5)+
    geom_hline(data = png, aes(yintercept=int), linetype="dashed", color = "lightgrey",size=0.5)+
    theme(strip.text.x = element_text(),
          strip.text.y = element_text(hjust = 0.5),
          strip.background = element_rect(color = 'black', linetype = 'solid'),
          strip.background.y = element_blank(),
          strip.background.x =element_blank(),
          panel.spacing=unit(1, "lines"),
          panel.background =element_rect(fill = 'white', colour = 'black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          legend.margin = margin(c(0.5, 2, 8, 25)),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black", size = 0, linetype = "solid"))
  
}

pdf(paste(plot_dir,'TFBS_effect_immune_cells.pdf',sep=''),width = 8,height = 4)
plot_allele_effect_tcells(strong_highfreq_effect)
dev.off()
