library(data.table);library(magrittr);library(dplyr)
library(reshape2);library(tidyr)
library(ComplexHeatmap)
library(wesanderson);library(RColorBrewer);library(circlize)
library(ggthemes);library(ggplot2);library(ggpubr)
library(R.utils)
library(openxlsx)

setwd('~/Desktop/Paper_1/D1_D2/')

## read snps 
read_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,chromatin_status := plyr::revalue(chrom_state,
                                           c('1_TssA'='active','2_TssAFlnk'='active','3_TxFlnk'='active','4_Tx'='active','5_TxWk'='active','6_EnhG'='active','7_Enh'='active','8_ZNF/Rpts'='active',
                                             '9_Het'='inactive','10_TssBiv'='inactive','11_BivFlnk'='inactive','12_EnhBiv'='inactive','13_ReprPC'='inactive','14_ReprPCWk'='inactive','15_Quies'='inactive'))
        ]
      
    )
  pop_names=c('denisova','neandertal','png') 
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

snps=read_snps('./')


read_snps=function(x){
  fread(x,sep=' ',header = T)[
    ,chromatin_status := plyr::revalue(chrom_state,
                                       c('1_TssA'='active','2_TssAFlnk'='active','3_TxFlnk'='active','4_Tx'='active','5_TxWk'='active','6_EnhG'='active','7_Enh'='active','8_ZNF/Rpts'='active',
                                         '9_Het'='inactive','10_TssBiv'='inactive','11_BivFlnk'='inactive','12_EnhBiv'='inactive','13_ReprPC'='inactive','14_ReprPCWk'='inactive','15_Quies'='inactive'))
    ]
}

denisova=read_snps('./denisovans_simplified.gz')


# enrichment(s)
enrichment=function(x,freq){
  if(freq=='nosplit'){
    # d0=copy(x)[deni_component%in%'0'] 
    d1=copy(x)[deni_component%in%'1'] 
    d2=copy(x)[deni_component%in%'2'] 
    # d1_d2=list(d1,d2)
    
    d1=d1[, .SD, .SDcols = !names(d1) %like% "freq"] %>% unique()
    d2=d2[, .SD, .SDcols = !names(d2) %like% "freq"]%>% unique()
    
      d1_vs_d2_enrichment=inner_join(d1,d2,by=c('chrom_state','chromatin_status','cell_type','cell_line')) %>%as.data.table()
d1_vs_d2_enrichment=d1_vs_d2_enrichment[
        ,enrichment :=fraction_snpdens_perstate_perepigenome.x/fraction_snpdens_perstate_perepigenome.y
        ][
          ,c('chrom_state','cell_type','cell_line','enrichment',
             'total_snps_perstate_perepigenome.x')
          ][,pop:='D1'] %>% unique()
    
    # pop_names=c('D1','D2') 
    # names(d1_d2)=pop_names
    # for (i in seq_along(d1_d2)){
    #   assign(pop_names[i],d1_d2[[i]],.GlobalEnv)}
    # d1_d2=Map(mutate,d1_d2,'pop'=names(d1_d2)) %>% lapply(function(z)setDT(z))
    
    return(d1_vs_d2_enrichment)
  } else{
    # d0=copy(x)[deni_component%in%'0'] 
    d1=copy(x)[deni_component%in%'1'] 
    d2=copy(x)[deni_component%in%'2'] 
    # d1_d2=list(d1,d2)
    
    d1=d1[, .SD, .SDcols = names(d1) %like% "freq|cell|chrom" ] %>% unique()
    d2=d2[, .SD, .SDcols = names(d2) %like% "freq|cell|chrom" ] %>% unique()
    
    d1_vs_d2_enrichment_freqsplit=inner_join(d1,d2,by=c('chrom_state','chromatin_status','cell_type','cell_line','freq_range')) %>%as.data.table()
d1_vs_d2_enrichment_freqsplit=d1_vs_d2_enrichment_freqsplit[
        ,enrichment :=fraction_snpdens_perstate_perepigenome_perfreqrange.x/fraction_snpdens_perstate_perepigenome_perfreqrange.y
        ][
          ,c('chrom_state','cell_type','cell_line','enrichment','freq_range',
             'total_snps_perstate_perepigenome_perfreqrange.x')
          ] [,pop:='D1'] %>% unique()

    # pop_names=c('D1','D2') 
    # names(d1_d2)=pop_names
    # for (i in seq_along(d1_d2)){
    #   assign(pop_names[i],d1_d2[[i]],.GlobalEnv)}
    # d1_d2=Map(mutate,d1_d2,'pop'=names(d1_d2)) %>% lapply(function(z)setDT(z))
    # 
    return(d1_vs_d2_enrichment_freqsplit)
  }
}

d1_vs_d2_nofreq=enrichment(denisova,'nosplit')
d1_vs_d2_freqsplit=enrichment(denisova,'split')

## matrix with log10 number of snps
numbsnps_matrix=function(x,split){
  make_matrix=function(x){
    df=copy(x) 
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
  if(split=='no'){
    df_nosplit=copy(x)
    df_nosplit=make_matrix(df_nosplit)
  }else if(split=='high'){
    df_high=copy(x) 
    df_high=df_high[freq_range%in%'high'][,freq_range:=NULL]
    df_high=make_matrix(df_high)
  } else{ 
    df_low=copy(x) 
    df_low=df_low[freq_range%in%'low'][,freq_range:=NULL]
    df_low=make_matrix(df_low)   
  }
}

numb_d1_nofreq=numbsnps_matrix(d1_vs_d2_nofreq,'no')
numb_d1_highfreq=numbsnps_matrix(d1_vs_d2_freqsplit,'high')
numb_d1_lowfreq=numbsnps_matrix(d1_vs_d2_freqsplit,'low')


# convert enrichment into matrix for heatmap
enrichment_matrix=function(x,split){
  build_matrix=function(x){
    df=copy(x) 
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
  
  if(split=='no'){
    df_nosplit=copy(x)
    df_nosplit=build_matrix(df_nosplit)
  }else if(split=='high'){
    df_high=copy(x) 
    df_high=df_high[freq_range%in%'high'][,freq_range:=NULL]
    df_high=build_matrix(df_high)
  } else{ 
    df_low=copy(x) 
    df_low=df_low[freq_range%in%'low'][,freq_range:=NULL]
    df_low=build_matrix(df_low)   
  }
}

d1_vs_d2_nofreq_matrix=enrichment_matrix(d1_vs_d2_nofreq,'no')
d1_vs_d2_highfreq_matrix=enrichment_matrix(d1_vs_d2_freqsplit,'high')
d1_vs_d2_lowfreq_matrix=enrichment_matrix(d1_vs_d2_freqsplit,'low')

## one sample t test for every cell ##
stat_test=function(x){
  df=copy(x)
  colnames(df)[5]='tot_numbsnps_perelement_perepigenome'
  df=df[,c('tot_numbsnps_perelement_perepigenome'):=NULL] %>% unique() 
  df=df[
    ,enrichment:=log2(enrichment)
    ]
  df=df %>% 
    split(as.factor(df$pop)) %>% 
    lapply(function(y)y %>% split(as.factor(y$chrom_state)) %>% 
             lapply(function(z) z %>%
                      mutate('test_stat' := t.test(z$enrichment,mu=0,alternative = 'g')$statistic,
                             'p_val':=t.test(z$enrichment,mu=0,alternative = 'g')$p.value,
                             'df' :=t.test(z$enrichment,mu=0,alternative = 'g')$parameter,
                             'p.adj':=p.adjust(p_val,method = 'bonferroni'),
                             'p.signif'= ifelse(`p.adj`<=0.0001,'****',
                                                ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                                       ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                                              ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))) %>% 
                      as.data.table()
             ) %>% rbindlist()
    )%>% rbindlist()
}



d1_vs_d2_nofreq_pvals=stat_test(d1_vs_d2_nofreq)
d1_vs_d2_highfreq_pvals=copy(d1_vs_d2_freqsplit)[freq_range%in%'high'][,freq_range:=NULL] %>% stat_test()
d1_vs_d2_lowfreq_pvals=copy(d1_vs_d2_freqsplit)[freq_range%in%'low'][,freq_range:=NULL] %>% stat_test()


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

pval_vector_nofreq=pval_vector(d1_vs_d2_nofreq_pvals) 
pval_vector_highfreq=pval_vector(d1_vs_d2_highfreq_pvals) 
pval_vector_lowfreq=pval_vector(d1_vs_d2_lowfreq_pvals) 


## make two Supp tables (1 enrichment and 1 pvals)
enrich_pval_tables=function(x,y,table){
  file_1=copy(x)
  file_1=file_1[order(as.factor(readr::parse_number(gsub("^.*\\.", "",file_1$chrom_state)))),]
  file_2=copy(y)
  file_2=file_2[order(as.factor(readr::parse_number(gsub("^.*\\.", "",file_2$chrom_state)))),]
  
  if(table=='pval'){
    pvalues=list(file_1,file_2)
    pvalues=lapply(pvalues,function(z)z=z[,c(2:4):=NULL] %>% unique())
    
    file_names=c('all_frequencies','common_to_high') 
    names(pvalues)=file_names
    for (i in seq_along(pvalues)){
      assign(file_names[i],pvalues[[i]],.GlobalEnv)}
    return(pvalues) 
  }else{
    
    enrichment=list(file_1,file_2)
    enrichment=lapply(enrichment,function(z)z=z[,c(1:5)] %>% unique())
    
    file_names=c('all_frequencies','common_to_high') 
    names(enrichment)=file_names
    for (i in seq_along(enrichment)){
      assign(file_names[i],enrichment[[i]],.GlobalEnv)}
    return(enrichment)
  }
  
}


enrichment_table=enrich_pval_tables(d1_vs_d2_nofreq_pvals,d1_vs_d2_highfreq_pvals,'enrichment')
write.xlsx(enrichment_table,'~/Desktop/Paper_1/pvalue_tables/Supplementary_tables/Supp_Table_5_D1vsD2_heatmap_enrichment.xlsx')

pval_table=enrich_pval_tables(d1_vs_d2_nofreq_pvals,d1_vs_d2_highfreq_pvals,'pval')
write.xlsx(pval_table,'~/Desktop/Paper_1/pvalue_tables/Supplementary_tables/Supp_Table_6_D1vsD2_heatmap_pvalues.xlsx')


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

chromatin_status=chromstatus_color(d1_vs_d2_nofreq_pvals)

chromatin_colors=chromatin_status$col
names(chromatin_colors)=chromatin_status$chrom_state

# pop
pop_color=function(x){
  x=x[,c('pop')
      ] %>% unique()
  x=x[
    ,col :=plyr::revalue(`pop`,c('D1'='burlywood1'))
    ]
}

pop_color=pop_color(d1_vs_d2_nofreq_pvals)
pop_colors=pop_color$col
names(pop_colors)=pop_color$pop

nihroadmap_colors=as.data.table(d1_vs_d2_nofreq_pvals)[
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


enrich_colors=colorRamp2(c(-2,-1,0,1,2), colors=c('blue3','royalblue3','white','red2','red3'))


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
                                                width = unit(5,'cm'),
                                                Cells = anno_simple(
                                                  gsub("^.*\\.", "",rownames(enrichmatrix)),
                                                  # height = unit(10,'cm'),
                                                  col= c("Adipose"='darkorange3',
                                                         "BCells"='palegreen4',
                                                         "Brain"='goldenrod3',
                                                         "Digestive"="plum3",
                                                         "Epithelial"="orange1",
                                                         "ES_cells"='hotpink4',
                                                         "ES_derived_cells"='dodgerblue2',
                                                         "Heart"="palevioletred2",
                                                         'IMR90_fetal_lung_fibroblast'='red3',
                                                         "iPSC"='mediumpurple4',
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

pdf('~/Desktop/Paper_1/heatmaps/d1_vs_d2_heatmap_enrichmet_nofreqsplit.pdf',width=25,height = 50)
enrichment_heatmap(d1_vs_d2_nofreq_matrix,numb_d1_nofreq,pval_vector_nofreq)
dev.off()


pdf('~/Desktop/Paper_1/heatmaps/d1_vs_d2_heatmap_enrichmet_highfreq.pdf',width=25,height = 50)
enrichment_heatmap(d1_vs_d2_highfreq_matrix,numb_d1_highfreq,pval_vector_highfreq)
dev.off()

# 
# pdf('~/Desktop/Paper_1/heatmaps/d1_d2_heatmap_enrichmet_lowfreq.pdf',width=25,height = 50)
# enrichment_heatmap(d1_d2_lowfreq_matrix,numb_d1_d2_lowfreq,pval_vector_lowfreq)
# dev.off()



# 
# lgd_chromHMM=Legend(title = '\n Chromatin states',
#                     labels=stringr::str_sort(names(chromHMM_colors),numeric = T),
#                     title_gp = gpar(fontsize=35,font=2),
#                     title_gap = unit(1, "cm"),
#                     title_position='topcenter',
#                     border = T,
#                     nrow = 3,
#                     grid_width = unit(1, "cm"),
#                     grid_height = unit(2, "cm"),
#                     gap=unit(0.5, "cm"),
#                     legend_gp = gpar(fill =c("#FF0000","#FF6E00","#32CD32","#008000","#006400","#C2E105", "#FFFF00","#66CDAA",
#                                              "#8A91D0","#CD5C5C","#E9967A","#BDB76B","#3A3838","#808080","#DCDCDC")),
#                    
#                     labels_gp = gpar(fontsize= 40))
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
# lgd_enrich=Legend(col_fun = enrich_colors, at = c(-2,-1,0,1,2),
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
# legends=packLegend(lgd_numbsnps,lgd_chromHMM,lgd_chromstatus,lgd_enrich,row_gap = unit(1, "cm"),direction = 'horizontal',gap = unit(10, "cm"))
# test=enrichment_heatmap(asnp_nofreq_matrix,numb_asnp_nofreq,pval_vector_nofreq)
# draw(legends)
# 
# pdf('~/Desktop/test.pdf',width=100,height = 200)
# draw(legends)
# dev.off()


