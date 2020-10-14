#### script to annotate SNPs with chromHMM data  
library(data.table);library(magrittr);library(dplyr);
library(GenomicRanges);library(R.utils)
library(wesanderson);library(RColorBrewer);
library(ggthemes);library(ggplot2);library(ggpubr)

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)
options(width=150)

### load cells 
setwd('/data/projects/punim0586/dvespasiani/Files/')

output_dir='./Archaic_introgression_in_PNG/Chromatin_states/SNPs_chromHMM_annotated/'
simplified_dir='./Archaic_introgression_in_PNG/Chromatin_states/simplified_set/'

plot_dir='./Archaic_introgression_in_PNG/Results/Plots/Chromatin_State/'

read_cells=function(x){
  cells=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y) 
      fread(y,sep = ' ',header=T) %>% 
        makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
        as.data.table()
    )
  cells=lapply(cells,function(x)x=x[,width:=as.numeric(width)][,state_coverage:=sum(width),by=.(chrom_state,cell_type)])
  cell_names=as.character(list.files(x,recursive=F,full.names=F))
  cell_names=gsub("\\..*","",cell_names)
  names(cells)=cell_names
  for (i in seq_along(cells)){
    assign(cell_names[i],cells[[i]],.GlobalEnv)}
  
  cells=cells[c(1:18)] # remove the encode ones
  
  return(cells)
  
}

cells=read_cells('./Annotation_and_other_files/Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states')

##  read snps
read_snps=function(x){
  pop=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)fread(y,sep=' ',header=T,
                            select = c('CHR','FROM','TO','REF','ALT',
                                       'ANC','all_frequency','freq_range')) %>% 
             setnames(c('seqnames','start','end','ref','alt','ancestral','all_freq','freq_range')) %>% 
             makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
             as.data.table() %>% unique()
    ) 
  pop_names=as.character(list.files(x,recursive=F,full.names=F))
  pop_names=gsub("\\..*","",pop_names)
  names(pop)=pop_names
  return(pop)
}

snps=read_snps('./Archaic_introgression_in_PNG/Grouped_filtered_snps/')

### merge snps with chromatin states ###
lapply(snps,function(x)setkey(x, seqnames, start, end))
lapply(cells,function(x)setkey(x, seqnames, start, end))


chromatin_state_annotation=function(cell,alleles){
  cell=lapply(cell,function(cell)
    alleles=lapply(alleles,function(alleles)
      foverlaps(cell,alleles, type="within")[
        ,element_seqnames:=seqnames
      ][
        ,c('seqnames','i.start','i.end','ref','alt','ancestral','all_freq','freq_range',
           'chrom_state','cell_type','state_coverage','element_seqnames','start','end')
          ] %>% setnames(old= c('start','end','i.start','i.end'),new = c('element_start','element_end','start','end'))
    )
  )
}

snps_chromatin_states=chromatin_state_annotation(snps,cells) 

denisova_chromatin_states=rbindlist(snps_chromatin_states[[1]])
neandertal_chromatin_states=rbindlist(snps_chromatin_states[[2]])
nonarchaics_chromatin_states=rbindlist(snps_chromatin_states[[3]])

cell_lines=function(x){x=x[,cell_line := plyr::revalue(cell_type,c("E017"="IMR90_fetal_lung_fibroblast", 
                                                                   
                                                                   "E002"="ES_cells","E008"="ES_cells","E001"="ES_cells",'E015'="ES_cells",'E014'="ES_cells",
                                                                   "E016"="ES_cells", "E003"='ES_cells',"E024"="ES_cells",
                                                                   
                                                                   "E020"="iPSC","E019"="iPSC","E018"="iPSC","E021"="iPSC","E022"="iPSC",
                                                                   
                                                                   "E007"="ES_derived_cells",
                                                                   "E009"="ES_derived_cells","E010"="ES_derived_cells","E013"="ES_derived_cells",
                                                                   "E012"="ES_derived_cells","E011"="ES_derived_cells","E004"="ES_derived_cells",
                                                                   "E005"="ES_derived_cells","E006"="ES_derived_cells",
                                                                   
                                                                   "E062"="TCells","E034"="TCells",
                                                                   "E045"="TCells","E033"="TCells","E044"="TCells","E043"="TCells",
                                                                   "E039"="TCells","E041"="TCells","E042"="TCells","E040"="TCells",
                                                                   "E037"="TCells","E048"="TCells","E038"="TCells","E047"="TCells",
                                                                   
                                                                   "E029"="BCells","E031"="BCells",
                                                                   "E035"="BCells","E051"="BCells","E050"="BCells","E036"="BCells",
                                                                   "E032"="BCells","E046"="BCells","E030"="BCells",
                                                                   
                                                                   "E026"="Mesenchymal","E049"="Mesenchymal",
                                                                   "E025"="Mesenchymal","E023"="Mesenchymal",
                                                                   
                                                                   "E052"='Myosatellite',
                                                                   
                                                                   "E055"="Epithelial","E056"="Epithelial","E059"="Epithelial",
                                                                   "E061"="Epithelial","E057"="Epithelial",
                                                                   "E058"="Epithelial","E028"="Epithelial","E027"="Epithelial",
                                                                   
                                                                   "E054"="Neurospheres","E053"="Neurospheres",
                                                                   
                                                                   "E112"="Thymus",'E093'="Thymus",
                                                                   
                                                                   "E071"="Brain","E074"="Brain",
                                                                   "E068"="Brain","E069"="Brain","E072"="Brain",
                                                                   "E067"="Brain","E073"="Brain","E070"="Brain",
                                                                   "E082"="Brain","E081"="Brain",
                                                                   
                                                                   "E063"="Adipose",
                                                                   
                                                                   "E100"="Muscle","E108"="Muscle","E107"="Muscle","E089"="Muscle","E090"="Muscle",
                                                                   
                                                                   "E083"="Heart","E104"="Heart","E095"="Heart","E105"="Heart","E065"="Heart",
                                                                   
                                                                   "E078"="Smooth_muscle","E076"="Smooth_muscle","E103"="Smooth_muscle","E111"="Smooth_muscle",
                                                                   
                                                                   "E092"="Digestive","E085"="Digestive","E084"="Digestive","E109"="Digestive",
                                                                   "E106"="Digestive","E075"="Digestive","E101"="Digestive","E102"="Digestive",
                                                                   "E110"="Digestive","E077"="Digestive","E079"="Digestive","E094"="Digestive",
                                                                   
                                                                   "E099"="Other_cells","E086"="Other_cells","E088"="Other_cells","E097"="Other_cells",
                                                                   "E087"="Other_cells","E080"="Other_cells",'E091'="Other_cells","E066"="Other_cells",
                                                                   "E098"="Other_cells", "E096"="Other_cells","E113"="Other_cells"
))]
}

snps_chromatin_states=list(denisova_chromatin_states,neandertal_chromatin_states,nonarchaics_chromatin_states) %>% 
  lapply(function(x)cell_lines(x))

names(snps_chromatin_states)=c('denisova','neandertal','png')

filenames_chromhmm=paste0(output_dir,names(snps_chromatin_states),sep='')
mapply(write.table,snps_chromatin_states, file = filenames_chromhmm,col.names = T, row.names = F, sep = " ", quote = F)

## count number snps per chromatin state 
asnps_enrichment=function(x,freq){
  df=copy(x)
  df=lapply(df,function(y)y=y[,numbsnps_perstate_percelltype:= .N, by=.(cell_type,chrom_state)
                              ][
                                ,c('chrom_state','cell_type','cell_line','numbsnps_perstate_percelltype')
                                ]%>% unique())
  
  enrichment=function(x,y){
  
    archaic=x[y,on=c('chrom_state','cell_type','cell_line')] %>% na.omit()
    
    archaic=archaic[
      ,asnps_nasnps_ratio:=numbsnps_perstate_percelltype/i.numbsnps_perstate_percelltype
      ][
        ,mean_asnps_nasnps_ratio:=mean(asnps_nasnps_ratio)
        ][
          ,mean_asnps_nasnps_ratio_chromstate:=mean(asnps_nasnps_ratio),by=.(chrom_state)
        ][
          ,enrichment:=asnps_nasnps_ratio/mean_asnps_nasnps_ratio
          ][
            ,c('chrom_state', 'cell_type','cell_line','enrichment','numbsnps_perstate_percelltype','asnps_nasnps_ratio','mean_asnps_nasnps_ratio','mean_asnps_nasnps_ratio_chromstate')
            ] %>% unique()
    
  }
  
  deni=copy(df[[1]])
  nean=copy(df[[2]])
  png=copy(df[[3]])
  
  deni_enrich=enrichment(deni,png)
  nean_enrich=enrichment(nean,png)
  archaics=list(deni_enrich,nean_enrich)
 
  if(freq=='high'){
    
    pop_names=c('denisova_high','neandertal_high')
    names(archaics)=pop_names
    for (i in seq_along( archaics)){
      assign(pop_names[i], archaics[[i]],.GlobalEnv)}
    return(archaics)
  }else{
    
   pop_names=c('denisova','neandertal')
  names( archaics)=pop_names
  for (i in seq_along( archaics)){
    assign(pop_names[i], archaics[[i]],.GlobalEnv)}
  return(archaics)
  }
  
}

snps_nofreqsplit=asnps_enrichment(snps_chromatin_states,'all')
snps_freqsplit=lapply(snps_chromatin_states,function(x)x=x[freq_range%in%'high'])
snps_freqsplit=asnps_enrichment(snps_freqsplit,'high')

## boxplot of asnps/nasnp ratio per chrom state 
ratio_plot=function(x){
  df=copy(x)
  df=lapply(df,function(x)x=x[,c(1:3,6,7)] %>% unique())
  df=Map(mutate,df,'pop'=names(df)) %>% rbindlist()
  df=df[order(as.factor(readr::parse_number(gsub("^.*\\.", "",df$chrom_state)))),][
    ,chrom_state:=gsub('.*_','',chrom_state)
    ]
  
  
  archaic_means=copy(df)[,c('mean_asnps_nasnps_ratio','pop')] %>% unique() 
  
  nihroadmap_colors=as.data.table(df)[
    ,col:=plyr::revalue(`chrom_state`,c('TssA'='#FF0000','TssAFlnk'='#FF6E00','TxFlnk'='#32CD32','Tx'='#008000',
                                        'TxWk'='#006400','EnhG'='#C2E105','Enh'='#FFFF00',
                                        'ZNF/Rpts'='#66CDAA','Het'='#8A91D0','TssBiv'='#CD5C5C','BivFlnk'='#E9967A',
                                        'EnhBiv'='#BDB76B','ReprPC'='#3A3838','ReprPCWk'='#808080',
                                        'Quies'='#DCDCDC'))
    ][
      ,c('chrom_state','col')
      ] %>% unique()
  

  my_palette_pop=c(#'#8D230F', #ambiguous
    '#C99E10', # denisova
    '#9B4F0F' #neandertal
  )
  
  names(my_palette_pop)= levels(as.factor(df$pop))
  colScale=scale_color_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal'))
  
  colState=scale_fill_manual(name= " ",values = nihroadmap_colors$col,labels = nihroadmap_colors$chrom_state)
  
  plot=ggplot(df,aes(x=chrom_state,y=asnps_nasnps_ratio,fill=chrom_state))+
    geom_violin(trim=F,scale = "width")+
  geom_boxplot(width=.4, position =  position_dodge(width = 0.4),outlier.size=0.2)+
     facet_wrap(pop~.,nrow=2)+
    geom_hline(data = archaic_means, aes(yintercept=mean_asnps_nasnps_ratio), linetype="dashed", color = my_palette_pop,size=1)+
    colScale+
    colState+
    ylab('\n aSNPs:naSNPs ratio per cell type')+
    xlab(' ')+
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
        legend.position = 'bottom',
        legend.margin = margin(c(0.5, 2, 8, 25)),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.y = element_text(hjust=0.5),
        axis.text=element_text(),
        axis.line = element_blank())
  return(plot)

}

##
pdf(paste(plot_dir,'snps_distribution/asnps_nasnps_ratio.pdf',sep=''),width=7,height =10)
ratio_plot(snps_nofreqsplit)
dev.off()

## barplot number of snps per chrom state
prepare_data=function(x){
  df=copy(x)
  df=lapply(df,function(x)x=x[,c(1:3,9)] %>% unique())
  df=lapply(df,function(x)x=x[
      ,tot_numbsnps_chromstate:=.N,by=.(chrom_state)
      ][
        ,log0_numbsnpchrom_state:=log10(tot_numbsnps_chromstate)
        ][,c('chrom_state','log0_numbsnpchrom_state')] %>% unique())
  
df=Map(mutate,df,'pop'=names(df)) %>% lapply(function(x)setDT(x))
}

numb_snps_chromstates=prepare_data(snps_chromatin_states)

log10_numbsnps_plot=function(df){
  x=copy(df)
  x=x[order(as.factor(readr::parse_number(gsub("^.*\\.", "",x$chrom_state)))),][
    ,chrom_state:=gsub('.*_','',chrom_state)
    ]
  nihroadmap_colors=as.data.table(x)[
    ,col:=plyr::revalue(`chrom_state`,c('TssA'='#FF0000','TssAFlnk'='#FF6E00','TxFlnk'='#32CD32','Tx'='#008000',
                                        'TxWk'='#006400','EnhG'='#C2E105','Enh'='#FFFF00',
                                        'ZNF/Rpts'='#66CDAA','Het'='#8A91D0','TssBiv'='#CD5C5C','BivFlnk'='#E9967A',
                                        'EnhBiv'='#BDB76B','ReprPC'='#3A3838','ReprPCWk'='#808080',
                                        'Quies'='#DCDCDC'))
    ][
      ,c('chrom_state','col')
      ] %>% unique()
  # colScale=scale_color_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal'))
  
  colState=scale_fill_manual(name= " ",values = nihroadmap_colors$col,labels = nihroadmap_colors$chrom_state)
  
  x$chrom_state=factor(x$chrom_state,levels = c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh', 'ZNF/Rpts',
                                                'Het','TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'))
  
    df_plot=ggplot(x,aes(x=chrom_state,y=log0_numbsnpchrom_state,fill=chrom_state))+
    geom_bar(stat='identity',position=position_dodge())+
      colState+
      xlab('')+
      ylab('\n Log10 number of aSNPs \n')+
    theme(
      legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      panel.background =element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(color = "black", size = 0.5, linetype = "solid"))
  
  return(df_plot)
  
}

pdf(paste(plot_dir,'snps_distribution/denisova_log10_numbsnps_chromstates.pdf',sep=''),width = 7,height =2)
log10_numbsnps_plot(numb_snps_chromstates[[1]])
dev.off()

pdf(paste(plot_dir,'snps_distribution/neandertal_log10_numbsnps_chromstates.pdf',sep=''),width = 7,height =2)
log10_numbsnps_plot(numb_snps_chromstates[[2]])
dev.off()

filenames_allfreq=paste0(paste0(simplified_dir,'all_freq/',sep=''),names(snps_nofreqsplit),sep='')
mapply(write.table,snps_nofreqsplit, file = filenames_allfreq,col.names = T, row.names = F, sep = " ", quote = F)

filenames_high=paste0(paste0(simplified_dir,'high_freq/',sep=''),names(snps_freqsplit),sep='')
mapply(write.table,snps_freqsplit, file = filenames_high,col.names = T, row.names = F, sep = " ", quote = F)

