## use this script to annotate the SNPs with chromatin state 
library(data.table)
library(magrittr)
library(dplyr)
library(GenomicRanges)
library(R.utils)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

options(width = 150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

input_dir = './filtered_files'
cells_dir = '../Annotation_and_other_files/Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states'
plot_dir = './Results/Plots/Chromatin_State/'
output_dir = './Chromatin_states/SNPs_chromHMM_annotated/'

table_dir='./Results/Tables/'

snps = list.files(input_dir,full.names = T,recursive = F,pattern='new_*') %>% lapply(function(y)fread(y,sep='\t',header = T))
names(snps)=c('denisova','modern_human','neanderthal')

## load cells
cells = list.files(cells_dir,recursive = F,full.names = T) %>% 
    lapply(function(y)fread(y,sep = ' ',header=T) %>% makeGRangesFromDataFrame(keep.extra.columns =T) %>% as.data.table()
)
cells = lapply(cells,function(x)x=x[,width:=as.numeric(width)][,state_coverage:=sum(width),by=.(chrom_state,cell_type)])
names(cells) = gsub("\\..*","",list.files(cells_dir,recursive=F,full.names=F))
cells = cells[c(1:18)] # remove the encode ones

### merge snps with chromatin states
lapply(snps,function(x)setkeyv(x,range_keys))
lapply(cells,function(x)setkeyv(x,range_keys))

chromatin_state_annotation = function(cell,alleles){
  cell=lapply(cell,function(cell)
    alleles=lapply(alleles,function(alleles)
      foverlaps(cell,alleles, type="within")[
        ,element_seqnames:=seqnames
      ][
        ,c('width','strand'):=NULL] %>%
        setnames(
          old= c(range_keys[-1],paste('i.',range_keys[-1],sep='')),
          new = c(paste('element',range_keys[-1],sep='_'),range_keys[-1])
          )
    )
  )
}

snps_chromatin_states = chromatin_state_annotation(snps,cells) 

denisova_chromatin_states = rbindlist(snps_chromatin_states[[1]])
human_chromatin_states = rbindlist(snps_chromatin_states[[2]])
neandertal_chromatin_states = rbindlist(snps_chromatin_states[[3]])

snps_chromatin_states = list(denisova_chromatin_states,human_chromatin_states,neandertal_chromatin_states) %>% lapply(function(x)group_celltypes(x))

names(snps_chromatin_states) = ancestry

# filenames_chromhmm = paste0(output_dir,'new_',names(snps_chromatin_states),'.txt.gz',sep='')
# mapply(fwrite,snps_chromatin_states, file = filenames_chromhmm,col.names = T, row.names = F, sep = "\t", quote = F)

## read this if steps above already done
# snps_chromatin_states = list.files(output_dir,full.names = T,recursive = F,pattern='new_*') %>% lapply(function(y)fread(y,sep='\t',header = T))
# names(snps_chromatin_states)=ancestry

count_snps_chromstate <- function(snps,select_cols){
  if('cell_line' %in% select_cols){
  counts <- copy(snps)
  counts<-counts[
      ,numbsnps_chromstate:=.N,by=.(chrom_state,cell_type)
      ][
        ,numbsnps:=.N,by=.(cell_type)
        ][
          ,numbsnps_notchromstate:=numbsnps-numbsnps_chromstate
          ]
          
  } else{
  counts <- copy(snps)
  counts<-counts[
      ,numbsnps_chromstate:=.N,by=.(chrom_state)
      ][
        ,numbsnps:=.N,
        ][
          ,numbsnps_notchromstate:=numbsnps-numbsnps_chromstate
          ]
         
  }
  final_counts <- copy(counts)%>%dplyr::select(c('numbsnps_chromstate','numbsnps_notchromstate',all_of(select_cols)))%>%unique()
  return(final_counts)
}

split_element_col <- function(x){
  x=copy(x)[,c('chrom_state','cell_type','cell_line'):= tstrsplit(elements, ".", fixed=TRUE)]
  return(x)
}

## genome wide enrichment (no cell type resolution)
snp_annot_counts <- copy(snps_chromatin_states)%>%lapply(function(x){
  x<-copy(x)[,c(..range_keys,'chrom_state')]%>%unique()%>%count_snps_chromstate(c('chrom_state'))
})

genwide_deni_or = calculate_or(snp_annot_counts[[1]],snp_annot_counts[[2]],c('chrom_state'))[,ancestry:='denisova']%>%adjust_pvalues()
genwide_nean_or = calculate_or(snp_annot_counts[[3]],snp_annot_counts[[2]],c('chrom_state'))[,ancestry:='neanderthal']%>%adjust_pvalues()

genwide_asnps_or <- rbind(genwide_deni_or,genwide_nean_or)

plot_or <- function(or){
  p<-ggplot(or, aes(x=factor(elements, level = chrom_states), y=odds_ratio,label = p_signif)) + 
    geom_point(aes(colour = elements))+
    geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
    scale_colour_manual(values = chromHMM_palette)+
    geom_hline(yintercept=1,linetype='dashed',size=.5)+
    geom_text(aes(y= or$upper_ci+0.1),size=5)+
    xlab(' ')+ylab('OR aSNPs vs naSNPs')+
    facet_wrap(~ancestry)+
    theme(
            panel.background =element_rect(fill = 'white', colour = 'black',size=1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.y = element_text(hjust = 0.5),
            strip.background = element_rect(color = 'black', linetype = 'solid'),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
            legend.position = "bottom",
            legend.key = element_rect(fill = "white", colour = "black"),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
            )
  return(p)
}

pdf(paste(plot_dir,'enrichment/genwide_asnps_or.pdf',sep=''),width=10,height = 7)
plot_or(genwide_asnps_or)
dev.off()


## plot mean numb snps per chrom state (mean across cell types)
meansnp_annot_counts <- copy(snps_chromatin_states[c(1,3)])%>%lapply(function(x){
  x<-copy(x)[
    ,c(..range_keys,'chrom_state','cell_line','cell_type')]%>%
    unique()%>%
    count_snps_chromstate(c('cell_line','cell_type','chrom_state'))
  x<-x[,numbsnps_notchromstate:=NULL][,meanNumbSnps:=round(mean(numbsnps_chromstate),1),by=.(chrom_state)
  ][
    ,c('chrom_state','meanNumbSnps')
    ]%>%unique()
})

meansnp_annot_counts <- Map(mutate,meansnp_annot_counts,ancestry=names(meansnp_annot_counts))%>%rbindlist()

# numbasnps_chromstate <- copy(snp_annot_counts[c(1,3)])
# numbasnps_chromstate <- Map(mutate,numbasnps_chromstate,ancestry=names(numbasnps_chromstate))%>%rbindlist()

plot_numbsnps <- function(x){
   p<-ggplot(x, aes(x=factor(chrom_state, level = chrom_states), y=log10(meanNumbSnps),fill=chrom_state)) + 
   geom_bar(stat='identity')+
    scale_fill_manual(values = chromHMM_palette)+
    xlab('')+ylab('log10 mean number aSNPs')+
    facet_wrap(~ancestry)+
    theme_classic()+
    theme(
            panel.background =element_rect(fill = 'white', colour = 'black',size=1),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.y = element_text(hjust = 0.5),
            strip.background.y = element_blank(),
            strip.background.x =element_blank(),
            legend.position = "none",
            legend.key = element_blank(),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
            )
  return(p)
}

pdf(paste(plot_dir,'enrichment/genwide_asnps_numbsnps_chromstate.pdf',sep=''),width=10,height = 3)
plot_numbsnps(meansnp_annot_counts)
dev.off()


highfreq_deni_or = calculate_or(highfreq_snp_annot_counts[[1]],highfreq_snp_annot_counts[[2]],c('chrom_state'))[,ancestry:='denisova']%>%adjust_pvalues()
highfreq_nean_or = calculate_or(highfreq_snp_annot_counts[[3]],highfreq_snp_annot_counts[[2]],c('chrom_state'))[,ancestry:='neanderthal']%>%adjust_pvalues()

highfreq_asnps_or <- rbind(highfreq_deni_or,highfreq_nean_or)

pdf(paste(plot_dir,'enrichment/highfreq_asnps_or.pdf',sep=''),width=10,height = 7)
plot_or(highfreq_asnps_or)
dev.off()

## repeat but for highfreq snps
mean_highfreq_snp_annot_counts <- copy(snps_chromatin_states[c(1,3)])%>%lapply(function(x){
  x<-copy(x)[!allfreq=='low'][
    ,c(..range_keys,'chrom_state','cell_line','cell_type')]%>%
    unique()%>%
    count_snps_chromstate(c('cell_line','cell_type','chrom_state'))
  x<-x[,numbsnps_notchromstate:=NULL][,meanNumbSnps:=round(mean(numbsnps_chromstate),1),by=.(chrom_state)
  ][
    ,c('chrom_state','meanNumbSnps')
    ]%>%unique()
})


mean_highfreq_snp_annot_counts <- Map(mutate,mean_highfreq_snp_annot_counts,ancestry=names(mean_highfreq_snp_annot_counts))%>%rbindlist()

pdf(paste(plot_dir,'enrichment/highfreq_asnps_numbsnps_chromstate.pdf',sep=''),width=10,height = 3)
plot_numbsnps(mean_highfreq_snp_annot_counts)
dev.off()

# ## enrichment chrom state/cell type combination

# ## now repeat the same steps but calculating the odds ratio at the cell type level
# snps_counts_celltypes <- copy(snps_chromatin_states)%>%lapply(function(x)
#   x<-copy(x)[,c(..range_keys,'chrom_state','cell_type','cell_line')]%>%unique()%>%count_snps_chromstate(c('chrom_state','cell_line','cell_type'))
# )

# denisovan_odds_ratio_celltype = calculate_or(snps_counts_celltypes[[1]],snps_counts_celltypes[[2]],c('chrom_state','cell_type','cell_line'))[,ancestry:='denisova']%>%adjust_pvalues()
# neanderthal_odds_ratio_celltype = calculate_or(snps_counts_celltypes[[3]],snps_counts_celltypes[[2]],c('chrom_state','cell_type','cell_line'))[,ancestry:='neanderthal']%>%adjust_pvalues()

# ##----------------------
# ## enrichment heatmap
# ##-----------------------
# ## Plotting both enrichment and depletion for this amount of data is tricky.
# ## best way is to set ns odds ratio pvals to 1 and then plot the (log2) of odds ratio 
# ## in this way the heatmap will have diverging colors starting from 0

# ## make matrix with the odd ratio enrichment value
# ## and with the number of SNPs annotated in each combination

# make_matrix = function(x,variable){
#   matrix = copy(x)[
#     ,celltype_line:=paste(cell_type,cell_line,sep='.')
#     ][
#       ,pop_chromstate:= paste(ancestry,chrom_state, sep='.')
#       ][
#         ,c('pop_chromstate','celltype_line','cell_type',..variable) 
#         ]%>%dcast(celltype_line ~ pop_chromstate, value.var=variable)
#   matrix[is.na(matrix)] = 1
#   rownames_matrix = matrix$celltype_line
#   matrix = matrix[,celltype_line:=NULL]%>%as.matrix()
#   rownames(matrix) = rownames_matrix
#   return(matrix)
# }

# test_enrichment = rbind(denisovan_odds_ratio_celltype,neanderthal_odds_ratio_celltype)[
#   ,c('chrom_state','cell_type','cell_line'):= tstrsplit(elements, ".", fixed=TRUE)
#   ][
#     ,odds_ratio:=ifelse(p<0.05,odds_ratio,1)
# ]%>%make_matrix('odds_ratio')%>%log2()
  
# test_numb_snps = copy(snps_counts_celltypes[c(1:2)])
# test_numb_snps = Map(mutate,test_numb_snps,ancestry = names(test_numb_snps))%>%rbindlist()%>%make_matrix('numbsnps_chromstate')%>%log10()

# ## color legends for heatmap
# # pop
# pop_color = my_palette_pop[c(1,3)]
# names(pop_color) = c('Denisova','Neanderthal')

# enrich_colors=colorRamp2(c(-3,-1,0,1,2), colors=c('blue3','royalblue1','white','red2','red3'))

# enrichment_heatmap=function(enrichmatrix,numbsnpsmatrix){
  
#   x = Heatmap(enrichmatrix, border = T,
#             rect_gp = gpar(type = "none"),
#             cell_fun = function(j, i, x, y, width, height, fill) {
#               grid.circle(x = x, y = y, r = abs(numbsnpsmatrix[i, j])/25, 
#                           gp = gpar(fill = enrich_colors(enrichmatrix[i, j]), col = 'black'))
              
#             },
#             row_dend_reorder = T,
#             row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
#             show_heatmap_legend = F,
#             show_row_names = F,
#             row_title =" ",
#             show_column_names = F,
#             column_order =order(as.factor(readr::parse_number(gsub("^.*\\.", "",colnames(enrichmatrix))))),
#             column_split = as.factor(readr::parse_number(gsub("^.*\\.", "", colnames(enrichmatrix)))),
#             column_title =' ',
#             top_annotation = HeatmapAnnotation(
#               which='col',
#               pop=anno_simple(sort(gsub("\\..*", "",colnames(enrichmatrix))),
#                               border = T,
#                               height=unit(2,'cm'),
#                               col=pop_color),
#               chrom_state=anno_simple(
#                 gsub("^.*\\.", "", colnames(enrichmatrix)),
#                 border=T,
#                 height = unit(2,'cm'),
#                 col=chromHMM_palette),
#               show_annotation_name = F),
#             left_annotation = HeatmapAnnotation(which='row',
#                                                 width = unit(6,'cm'),
#                                                 Cells = anno_simple(
#                                                   gsub("^.*\\.", "",rownames(enrichmatrix)),
#                                                   col= tissue_colors),
#                                                   show_annotation_name = F))
#   return(x)
# }

# pdf(paste(plot_dir,'enrichment/heatmap_all_asnp_oddsratio.pdf',sep=''),width=25,height = 50)
# enrichment_heatmap(test_enrichment,test_numb_snps)
# dev.off()

# ## get legends
# lgd_numbsnps=Legend(
#   at=c('1','2','3','4','5'),
#   title = '\n log10 number aSNPs',
#   title_position='topcenter',
#   title_gp = gpar(fontsize=35,font=2),
#   title_gap = unit(1,'cm'),
#   type='points',
#   pch=1,
#   nrow=5,
#   background = "white",
#   legend_gp = gpar(fill ='white',fontsize=35),
#   grid_width = unit(0.8, "cm"),
#   grid_height = unit(2, "cm"),
#   labels_gp = gpar(fontsize= 40)
# )
 
# lgd_enrich=Legend(
#   col_fun = enrich_colors, at = c(-2,-1,0,1,2),
#   title='\n log2 aSNPs OR',
#   direction='horizontal',
#   border=T,
#   title_position='topcenter',
#   title_gp = gpar(fontsize=35,font=2),
#   title_gap = unit(1,'cm'),
#   legend_height = unit(10, "cm"),
#   grid_width = unit(0.8, "cm"),
#   grid_height = unit(2, "cm"),
#   labels_gp = gpar(fontsize= 40),
#   gap = unit(4, "cm")
# )

# lgd_chromHMM=Legend(
#   title = '\n Chromatin states',
#   labels=stringr::str_sort(names(chromHMM_palette),numeric = T),
#   title_gp = gpar(fontsize=35,font=2),
#   title_gap = unit(1, "cm"),
#   title_position='topcenter',
#   border = T,
#   nrow = 3,
#   grid_width = unit(1, "cm"),
#   grid_height = unit(2, "cm"),
#   gap=unit(0.5, "cm"),
#   legend_gp = gpar(
#     fill = chromHMM_palette
#       ),
#   labels_gp = gpar(fontsize= 40)
# )

# legends=packLegend(
#   lgd_enrich,lgd_chromHMM,lgd_numbsnps,
#   row_gap = unit(1, "cm"),
#   direction = 'horizontal',
#   gap = unit(10, "cm")
# ) 

# pdf(paste(plot_dir,'heatmap_allsnps_legends.pdf',sep=''),width=100,height = 200)
# draw(legends)
# dev.off()


## write spreadsheet with OR results
## 1) genome-wide enrichment allfreq and highfreq
## 2) enrichment immune cells chromstate allfreq 
## 3) enrichment immune cells chromstate highfreq

numbsnps_genwide <- lapply(snps_chromatin_states,function(x){
  x<-copy(x)[,c(..range_keys,'chrom_state')]%>%unique()
  x<-x[,numbsnps:=.N,by=.(chrom_state)][,c('numbsnps','chrom_state')]%>%unique()
  }
)
numb_asnps_genwide <- copy(numbsnps_genwide[c(1,3)])
numb_asnps_genwide <- Map(mutate,numb_asnps_genwide,ancestry=names(numb_asnps_genwide))%>%rbindlist()

supp_table_or_genomewide <- copy(genwide_asnps_or)%>%setnames(old='elements',new='chrom_state')
supp_table_or_genomewide <- supp_table_or_genomewide[numb_asnps_genwide,on=c('chrom_state','ancestry'),nomatch=0]
supp_table_or_genomewide <- supp_table_or_genomewide[order(as.factor(readr::parse_number(supp_table_or_genomewide$chrom_state))),]

## highfreq snps
numb_highfreq_asnps <- copy(highfreq_snp_annot_counts[c(1,3)])
numb_highfreq_asnps <- Map(mutate,numb_highfreq_asnps,ancestry=names(numb_highfreq_asnps))%>%rbindlist()
numb_highfreq_asnps <- numb_highfreq_asnps[,numbsnps_notchromstate:=NULL]%>%setnames(old='numbsnps_chromstate',new='numb_snps')

supp_table_or_highfreq <- copy(highfreq_asnps_or)%>%setnames(old='elements',new='chrom_state')
supp_table_or_highfreq <- supp_table_or_highfreq[numb_highfreq_asnps,on=c('ancestry','chrom_state'),nomatch=0]
supp_table_or_highfreq <- supp_table_or_highfreq[order(as.factor(readr::parse_number(supp_table_or_highfreq$chrom_state))),]


spreadsheet = list(supp_table_or_genomewide,supp_table_or_highfreq)
names(spreadsheet) = c('genome wide chrom state OR','commtohigh freq chrom state OR')

write.xlsx(spreadsheet,paste(table_dir,'Supp_Table_OR_chromstate_pvals.xlsx',sep=''),append=T,overwrite=T)

