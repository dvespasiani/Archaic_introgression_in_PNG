## This script looks at the pleiotropic activiity of the functional elements carrying SNPs for each chrom state

library(data.table)
library(magrittr)
library(dplyr)
library(GenomicRanges)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(R.utils)
library(openxlsx)

options(width = 150)


### load cells 
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

input_dir = './Chromatin_states/SNPs_chromHMM_annotated'
# table_dir='./Results/Tables/'
plot_dir='./Results/Plots/Chromatin_State/'

range_keys = c('seqnames','start','end')
columns_to_read = c(range_keys,'REF','ALT','chrom_state','cell_type','cell_line')

snps_chrom_states = list.files(input_dir,recursive = F,full.names = T) %>%
  lapply(function(y)fread(y,sep='\t',header=T,select = columns_to_read)%>% unique()
)
names(snps_chrom_states) = gsub("\\..*","",list.files(input_dir,recursive = F,full.names = F))


## calculate pleiotropy
## pleiotropy here is intendeted out of 111 cell types in how many cell types each SNP has same functional annotation
chrom_state_pleiotropy = copy(snps_chrom_states)%>%
  lapply(function(x)x=x[
    ,pleiotropy:=.N,by=.(seqnames,start,end,chrom_state)
    ][
      ,numb_snps_chromstate := .N ,by=.(chrom_state)
    ][
      ,numb_snps_chromstate_pleiotropy  := .N ,by=.(pleiotropy,chrom_state)
    ][
      ,propsnps_chromstate_pleiotropy := numb_snps_chromstate_pleiotropy/numb_snps_chromstate
    ][
    ,c('pleiotropy','propsnps_chromstate_pleiotropy','chrom_state')
    ]%>%unique()%>%setorderv('pleiotropy',1)
)

chrom_state_pleiotropy = lapply(
  chrom_state_pleiotropy,
  function(x)x=x%>%split(by='chrom_state')%>%
  lapply(
    function(y)y=y[,cumsum := cumsum(propsnps_chromstate_pleiotropy)]
    )%>%rbindlist()
)

chrom_state_levels = c(
      '1_TssA','2_TssAFlnk','3_TxFlnk','4_Tx','5_TxWk','6_EnhG','7_Enh', '8_ZNF/Rpts',
      '9_Het','10_TssBiv','11_BivFlnk','12_EnhBiv','13_ReprPC','14_ReprPCWk','15_Quies'
)

chrom_state_pleiotropy = Map(mutate,chrom_state_pleiotropy,pop=names(chrom_state_pleiotropy))%>%
  rbindlist()%>%
  arrange(
    factor(chrom_state, levels = chrom_state_levels)
)


## plot the cumulative distribution of the proportion of SNPs across the epigenomes

my_palette_pop=c(
    '#C99E10', # denisova
    '#1E656D', #png
    '#9B4F0F' #neandertal
)

pleiotropy_plot = function(df){

  df$chrom_state = factor(df$chrom_state,levels = chrom_state_levels)
  
  plot = ggplot(df,aes(x = pleiotropy,y = cumsum,color = pop))+
  geom_line()+
  scale_color_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Modern Humans','Neanderthal'))+
  facet_wrap(chrom_state~., ncol = 3)+
  xlab('Number of epigenomes')+ylab('\n Cumulative proportion SNPs \n ')+
  theme(
    strip.text.x = element_text(),
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
    axis.line = element_blank()
    )
  return(plot)

}


pdf(paste(plot_dir,'pleiotropy/cumulative_pleiotropy.pdf',sep=''),width = 7,height = 7)
pleiotropy_plot(chrom_state_pleiotropy)
dev.off()
