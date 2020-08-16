library(data.table);library(magrittr);library(dplyr)
library(reshape2);library(tidyr)
library(ComplexHeatmap)
library(wesanderson);library(RColorBrewer);library(circlize)
library(ggthemes);library(ggplot2);library(ggpubr)
library(R.utils)
library(openxlsx)

setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG')

simplified_dir='./Chromatin_states/simplified_set/new_set/'

plot_dir='./Results/Plots/Chromatin_State/'
table_dir='./Results/Tables/'

## read snps 
read_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,chromatin_status := plyr::revalue(chrom_state,
                                           c('1_TssA'='active','2_TssAFlnk'='active','3_TxFlnk'='active','4_Tx'='active','5_TxWk'='active','6_EnhG'='active','7_Enh'='active','8_ZNF/Rpts'='active',
                                             '9_Het'='inactive','10_TssBiv'='inactive','11_BivFlnk'='inactive','12_EnhBiv'='inactive','13_ReprPC'='inactive','14_ReprPCWk'='inactive','15_Quies'='inactive'))
        ][
          ,cell_line:=ifelse(cell_line=='iPSC','IPSC',cell_line)
        ][
          ,cell_line:=ifelse(cell_line=='ES_derived_cells','Es_derived_cells',cell_line)
          ][
            ,cell_line:=ifelse(cell_line=='ES_cells','Es_cells',cell_line)
            ]
      
    )
  pop_names=c('denisova','neandertal') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  x=Map(mutate,x,'pop'=names(x)) %>% lapply(function(z)setDT(z))
  
  return(x)
}

asnp_nofreq=read_snps(paste(simplified_dir,'all_freq/',sep=''))
asnp_freq=read_snps(paste(simplified_dir,'high_freq/',sep=''))

## matrix with log10 number of snps
numbsnps_matrix=function(x){
  df=copy(x) %>% rbindlist()
  colnames(df)[5]='tot_numbsnps_perelement_perepigenome'
  df=df[
    ,c('chrom_state','cell_type','cell_line','tot_numbsnps_perelement_perepigenome','pop')
    ] %>% unique()
  df=df[, cell_type_line := paste(cell_type,cell_line, sep='.')
        ][
          ,pop_chromstate := paste(pop,chrom_state, sep='.')
          ][
            ,c('cell_type','cell_line','pop','chrom_state'):=NULL
            ]%>% unique() 
  df=dcast(df,cell_type_line ~  pop_chromstate,value.var='tot_numbsnps_perelement_perepigenome')
  rownames(df)=df$cell_type
  df=df%>% dplyr::select(-1) 
  df=df%>% as.matrix() 
  df=df%>% log10()
  
}

numb_asnp_nofreq=numbsnps_matrix(asnp_nofreq)
numb_asnp_highfreq=numbsnps_matrix(asnp_freq)

# convert enrichment into matrix for heatmap
enrichment_matrix=function(x){
  df=copy(x) %>% rbindlist()
  colnames(df)[5]='tot_numbsnps_perelement_perepigenome'
  
  df=df[, cell_type_line := paste(cell_type,cell_line, sep='.')
        ][
          ,pop_chromstate := paste(pop,chrom_state, sep='.')
          ][
            ,c('cell_type_line','pop_chromstate','enrichment')
            ]%>% unique() 
  df=dcast(df,cell_type_line ~  pop_chromstate,value.var='enrichment')
  rownames(df)=df$cell_type
  df=df%>% dplyr::select(-1)
  df=df%>% as.matrix() %>% log2()
}


asnp_nofreq_matrix=enrichment_matrix(asnp_nofreq)
asnp_highfreq_matrix=enrichment_matrix(asnp_freq)

## one sample t test for every cell ##
stat_test=function(x){
  df=copy(x)
  colnames(df)[5]='tot_numbsnps_perelement_perepigenome'
  df=df[,c('tot_numbsnps_perelement_perepigenome'):=NULL] %>% unique() 
  df=df %>% 
    split(as.factor(df$pop)) %>% 
    lapply(function(y)y %>% split(as.factor(y$chrom_state)) %>% 
             lapply(function(z) z %>%
                      mutate('test_stat' := t.test(z$enrichment,mu= 1,alternative = 'g')$statistic,
                             'p':=t.test(z$enrichment,mu= 1,alternative = 'g')$p.value,
                             'df' :=t.test(z$enrichment,mu= 1,alternative = 'g')$parameter)%>% 
                      as.data.table()
             ) %>% rbindlist()
    )%>% rbindlist()
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


asnp_nofreq_pvals=lapply(asnp_nofreq,function(x)stat_test(x) %>% adjust_pvalues()) %>% rbindlist()
asnp_highfreq_pvals=lapply(asnp_freq,function(x)stat_test(x)%>% adjust_pvalues())%>% rbindlist()


pval_vector=function(x){
  x=x[,c('chrom_state','p.adj','pop','test_stat')
      ][
        ,pop_chromstate := paste(pop,chrom_state, sep='.')
        ] %>% unique()
  x=x[,c('chrom_state','pop'):=NULL
      ][
        ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                           ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                  ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                         ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
        ] %>% arrange(`pop_chromstate`)
  
  named_vector=structure(as.character(x$p.signif),names=as.character(x$pop_chromstate)) 
  
}

pval_vector_nofreq=pval_vector(asnp_nofreq_pvals) 
pval_vector_highfreq=pval_vector(asnp_highfreq_pvals)

## make two Supp tables (1 enrichment and 1 pvals)
enrich_pval_tables=function(x,y,table){
  file_1=copy(x)
  file_1=file_1[order(as.factor(readr::parse_number(gsub("^.*\\.", "",file_1$chrom_state)))),]
  file_2=copy(y)
  file_2=file_2[order(as.factor(readr::parse_number(gsub("^.*\\.", "",file_2$chrom_state)))),]
  
  if(table=='pval'){
    pvalues=list(file_1,file_2)
    pvalues=lapply(pvalues,function(z)z=z[,c(2:6):=NULL] %>% unique())
    
    file_names=c('all_frequencies','common_to_high') 
    names(pvalues)=file_names
    for (i in seq_along(pvalues)){
      assign(file_names[i],pvalues[[i]],.GlobalEnv)}
    return(pvalues) 
  }else{
    
    enrichment=list(file_1,file_2)
    enrichment=lapply(enrichment,function(z)z=z[,c(1:5,7,9)] %>% unique())
    
    file_names=c('all-frequencies','common-to-high') 
    names(enrichment)=file_names
    for (i in seq_along(enrichment)){
      assign(file_names[i],enrichment[[i]],.GlobalEnv)}
    return(enrichment)
  }
  
}


enrichment_table=enrich_pval_tables(asnp_nofreq_pvals,asnp_highfreq_pvals,'enrichment')
write.xlsx(enrichment_table,paste(table_dir,'Supp_Table_2_heatmap_enrichment.xlsx',sep=''))

pval_table=enrich_pval_tables(asnp_nofreq_pvals,asnp_highfreq_pvals,'pval')
write.xlsx(pval_table,paste(table_dir,'Supp_Table_3_heatmap_pvalues.xlsx',sep=''))

## color legends for heatmap
chromstatus_color=function(x){
  x=x[
    ,c('chrom_state')
    ] %>% unique()
  x=x[
    ,col :=plyr::revalue(`chrom_state`,c('1_TssA'='khaki3','2_TssAFlnk'='khaki3','3_TxFlnk'='khaki3',
                                         '4_Tx'='khaki3','5_TxWk'='khaki3','6_EnhG'='khaki3',
                                         '7_Enh'='khaki3','8_ZNF/Rpts'='khaki3','9_Het'='ivory3',
                                         '10_TssBiv'='ivory3','11_BivFlnk'='ivory3','12_EnhBiv'='ivory3',
                                         '13_ReprPC'='ivory3','14_ReprPCWk'='ivory3','15_Quies'='ivory3'))
    ]
}

chromatin_status=chromstatus_color(asnp_nofreq_pvals)

chromatin_colors=chromatin_status$col
names(chromatin_colors)=chromatin_status$chrom_state

# pop
pop_color=function(x){
  x=x[,c('pop')
      ] %>% unique()
  x=x[   
    ,col :=plyr::revalue(`pop`,c('denisova'='#C99E10','neandertal'='#9B4F0F'))
    ]
}

pop_color=pop_color(asnp_nofreq_pvals)
pop_colors=pop_color$col
names(pop_colors)=pop_color$pop

nihroadmap_colors=as.data.table(asnp_nofreq_pvals)[
  ,col:=plyr::revalue(`chrom_state`,c('1_TssA'='#FF0000','2_TssAFlnk'='#FF6E00','3_TxFlnk'='#32CD32','4_Tx'='#008000',
                                      '5_TxWk'='#006400','6_EnhG'='#C2E105','7_Enh'='#FFFF00',
                                      '8_ZNF/Rpts'='#66CDAA','9_Het'='#8A91D0','10_TssBiv'='#CD5C5C','11_BivFlnk'='#E9967A',
                                      '12_EnhBiv'='#BDB76B','13_ReprPC'='#3A3838','14_ReprPCWk'='#808080',
                                      '15_Quies'='#DCDCDC'))
  ][
    ,c('chrom_state','col')
    ] %>% unique()

chromHMM_colors=nihroadmap_colors$col
names(chromHMM_colors)=nihroadmap_colors$chrom_state


enrich_colors=colorRamp2(c(-5,-2,-1,0,1,2), colors=c('blue3','royalblue3','royalblue1','white','red2','red3'))


enrichment_heatmap=function(enrichmatrix,numbsnpsmatrix,pval){
  
  x=Heatmap(enrichmatrix, border = T,
            rect_gp = gpar(type = "none"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.circle(x = x, y = y, r = abs(numbsnpsmatrix[i, j])/25,  #asnp_chromstates[i, j])/5 * min(unit.c(width, height)
                          gp = gpar(fill = enrich_colors(enrichmatrix[i, j]), col = 'black'))
              
            },
            # height=unit(20,'cm'),
            row_dend_reorder = T,
            row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
            #row_split = as.factor(gsub("^.*\\.", "", rownames(asnpenrichment_nofreqsplit))),
            show_heatmap_legend = F,
            show_row_names = F,
            row_title =" ",
            show_column_names = F,
            column_order =order(as.factor(readr::parse_number(gsub("^.*\\.", "",colnames(enrichmatrix))))),
            column_split = as.factor(readr::parse_number(gsub("^.*\\.", "", colnames(enrichmatrix)))),
            column_title =' ',
            top_annotation = HeatmapAnnotation(
              which='col',
              pop=anno_simple(sort(gsub("\\..*", "",colnames(enrichmatrix))),
                              border = T,
                              height=unit(2,'cm'),
                              col=pop_colors),
              new_chromatin_status=anno_simple(
                gsub("^.*\\.", "",colnames(enrichmatrix)),
                border = T,
                height = unit(2,'cm'),
                col=chromatin_colors,
                pch = pval,
                pt_gp=gpar(fontsize=40)),
              chrom_state=anno_simple(
                gsub("^.*\\.", "", colnames(enrichmatrix)),
                border=T,
                height = unit(2,'cm'),
                col=chromHMM_colors),
              show_annotation_name = F),
            left_annotation = HeatmapAnnotation(which='row',
                                                width = unit(6,'cm'),
                                                Cells = anno_simple(
                                                  gsub("^.*\\.", "",rownames(enrichmatrix)),
                                                  # height = unit(10,'cm'),
                                                  col= c("Adipose"='darkorange3',
                                                         "BCells"='palegreen4',
                                                         "Brain"='goldenrod3',
                                                         "Digestive"="plum3",
                                                         "Epithelial"="orange1",
                                                         "Es_cells"='hotpink4',
                                                         "Es_derived_cells"='dodgerblue2',
                                                         "Heart"="palevioletred2",
                                                         'IMR90_fetal_lung_fibroblast'='red3',
                                                         "IPSC"='mediumpurple4',
                                                         "Mesenchymal"='hotpink3',
                                                         "Muscle"="indianred",
                                                         "Myosatellite"='darkorange2',
                                                         "Neurospheres"='lightgoldenrod2',
                                                         "Other_cells"="grey60",
                                                         "Smooth_muscle"="hotpink1",
                                                         "TCells"='chartreuse4',
                                                         "Thymus"='yellow2')),show_annotation_name = F))
  return(x)
}

pdf(paste(plot_dir,'enrichment/heatmap_asnpenrichmet_nofreqsplit.pdf',sep=''),width=25,height = 50)
enrichment_heatmap(asnp_nofreq_matrix,numb_asnp_nofreq,pval_vector_nofreq)
dev.off()


pdf(paste(plot_dir,'enrichment/heatmap_asnpenrichmet_highfreq.pdf',sep=''),width=25,height = 50)
enrichment_heatmap(asnp_highfreq_matrix,numb_asnp_highfreq,pval_vector_highfreq)
dev.off()

# 
# lgd_chromstatus=Legend(title = '\n State activity ', labels = c("Active","Inactive"),
#                        title_gp = gpar(fontsize=35,font=2),
#                        title_gap = unit(1,'cm'),
#                        border = T,nrow = 2,
#                        gap=unit(3, "cm"),
#                        grid_width = unit(0.8, "cm"),
#                        grid_height = unit(2, "cm"),
#                        legend_gp = gpar(fill = c("Active"='khaki3', "Inactive"='ivory3')),
#                        labels_gp = gpar(fontsize= 40))
# 
# 
# lgd_numbsnps=Legend(at=c('1','2','3','4','5'),
#                     title = '\n log10 number aSNPs',
#                     title_position='topcenter',
#                     title_gp = gpar(fontsize=35,font=2),
#                     title_gap = unit(1,'cm'),
#                     type='points',
#                     pch=1,
#                     nrow=5,
#                     background = "white",
#                     legend_gp = gpar(fill ='white',fontsize=35),
#                     grid_width = unit(0.8, "cm"),
#                     grid_height = unit(2, "cm"),
#                     labels_gp = gpar(fontsize= 40))
# 
# 
# enrich_colors_lgd=colorRamp2(c(-2,-1,0,1,2), colors=c('blue2','royalblue3','white','red','red1'))
# 
# lgd_enrich=Legend(col_fun = enrich_colors_lgd, at = c(-2,-1,0,1,2),
#                   title='\n log2 aSNPs enrichment',
#                   direction='horizontal',
#                   border=T,
#                   title_position='topcenter',
#                   title_gp = gpar(fontsize=35,font=2),
#                   title_gap = unit(1,'cm'),
#                   legend_height = unit(10, "cm"),
#                   grid_width = unit(0.8, "cm"),
#                   grid_height = unit(2, "cm"),
#                   labels_gp = gpar(fontsize= 40),
#                   gap = unit(4, "cm"))
# 
# legends=packLegend(lgd_enrich,row_gap = unit(1, "cm"),direction = 'horizontal',gap = unit(10, "cm"))
# grid.draw(legends)
# 
# pdf('~/Desktop/test.pdf',width=100,height = 200)
# draw(legends)
# dev.off()


