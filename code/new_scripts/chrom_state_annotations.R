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

snps = list.files(input_dir,full.names = T,recursive = F) %>% lapply(function(y)fread(y,sep='\t',header = T))
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

snps_chromatin_states = list(denisova_chromatin_states,neandertal_chromatin_states,human_chromatin_states) %>% lapply(function(x)group_celltypes(x))

names(snps_chromatin_states) = c('Denisova','Neanderthal','Modern_humans')

filenames_chromhmm = paste0(output_dir,names(snps_chromatin_states),'.txt',sep='')
mapply(write.table,snps_chromatin_states, file = filenames_chromhmm,col.names = T, row.names = F, sep = "\t", quote = F)

## compute odds ratio aSNPs per chromatin state
snps_annotation_counts = copy(snps_chromatin_states)%>%lapply(function(x)x%>%dplyr::select(-c(contains('element'),'state_coverage',contains('cell')))%>%unique())
snps_annotation_counts = lapply(snps_annotation_counts,function(x)x=x[,numb_snsp_chromstate := .N ,by=.(chrom_state)][,numb_snps:=.N][,numb_snps_notin_chromstate:=numb_snps-numb_snsp_chromstate][,c('chrom_state','numb_snsp_chromstate','numb_snps_notin_chromstate')]%>%unique())

denisovan_odds_ratio = calculate_odds_ratio(snps_annotation_counts[[1]],snps_annotation_counts[[3]],'chrom_state',T)[,ancestry:='Denisova']%>%adjust_pvalues()
neanderthal_odds_ratio = calculate_odds_ratio(snps_annotation_counts[[2]],snps_annotation_counts[[3]],'chrom_state',T)[,ancestry:='Neanderthal']%>%adjust_pvalues()

asnps_odds_ratio = rbind(denisovan_odds_ratio,neanderthal_odds_ratio)
asnps_odds_ratio = asnps_odds_ratio[
  ,p_signif:=ifelse(p_signif==' ',' ','*')
]
## plot results
chromstate_levels = unique(asnps_odds_ratio$elements)%>%stringr::str_sort( numeric = TRUE)
chromHMM_palette = c(
  '#FF0000','#FF6E00','#32CD32',
  '#008000','#006400','#C2E105',
  '#FFFF00','#66CDAA','#8A91D0',
  '#CD5C5C','#E9967A','#BDB76B',
  '#3A3838','#808080','#DCDCDC'
)
names(chromHMM_palette) = chromstate_levels


pdf(paste(plot_dir,'chromstate_odds_ratio_asnps_vs_nasnps.pdf',sep=''),width=8,height = 5)
ggplot(asnps_odds_ratio, aes(x=factor(elements, level = chromstate_levels), y=odds_ratio,label = p_signif)) + 
geom_point(aes(colour = elements))+
geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
scale_colour_manual(values = chromHMM_palette)+
geom_hline(yintercept=1,linetype='dashed',size=.5)+
geom_text(aes(y= asnps_odds_ratio$upper_ci+0.1),size=5)+
xlab(' ')+ylab('odds ratio \n aSNPs vs naSNPs')+
facet_wrap(~ancestry)+
theme(
        panel.background =element_rect(fill = 'white', colour = 'black',size=1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(hjust = 0.5),
        strip.background = element_rect(color = 'black', linetype = 'solid'),
        strip.background.y = element_blank(),
        strip.background.x =element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white", colour = "black"),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()


## now repeat the same steps but calculating the odds ratio at the cell type level
snps_counts_celltypes = copy(snps_chromatin_states)%>%lapply(function(x)x%>%dplyr::select(-c(contains('element'),'state_coverage'))%>%unique())
snps_counts_celltypes = lapply(
  snps_counts_celltypes,function(x)
  x=x[
    ,numb_snsp_chromstate_celltype := .N ,by=.(chrom_state,cell_type)
    ][
      ,numb_snps_celltype:=.N,by=.(cell_type)
      ][
        ,numb_snps_chromstate_notin_celltype:=numb_snps_celltype-numb_snsp_chromstate_celltype
        ][
          ,c('chrom_state','cell_type','cell_line','numb_snsp_chromstate_celltype','numb_snps_chromstate_notin_celltype')
          ]%>%unique()
)


denisovan_odds_ratio_celltype = calculate_odds_ratio(snps_counts_celltypes[[1]],snps_counts_celltypes[[3]],c('chrom_state','cell_type','cell_line'),T)[,ancestry:='Denisova']%>%adjust_pvalues()
neanderthal_odds_ratio_celltype = calculate_odds_ratio(snps_counts_celltypes[[2]],snps_counts_celltypes[[3]],c('chrom_state','cell_type','cell_line'),T)[,ancestry:='Neanderthal']%>%adjust_pvalues()

##----------------------
## enrichment heatmap
##-----------------------
## Plotting both enrichment and depletion for this amount of data is tricky.
## best way is to set ns odds ratio pvals to 1 and then plot the (log2) of odds ratio 
## in this way the heatmap will have diverging colors starting from 0

## make matrix with the odd ratio enrichment value
## and with the number of SNPs annotated in each combination

make_matrix = function(x,variable){
  matrix = copy(x)[
    ,celltype_line:=paste(cell_type,cell_line,sep='.')
    ][
      ,pop_chromstate:= paste(ancestry,chrom_state, sep='.')
      ][
        ,c('pop_chromstate','celltype_line','cell_type',..variable) 
        ]%>%dcast(celltype_line ~ pop_chromstate, value.var=variable)
  matrix[is.na(matrix)] = 1
  rownames_matrix = matrix$celltype_line
  matrix = matrix[,celltype_line:=NULL]%>%as.matrix()
  rownames(matrix) = rownames_matrix
  return(matrix)
}

test_enrichment = rbind(denisovan_odds_ratio_celltype,neanderthal_odds_ratio_celltype)[
  ,c('chrom_state','cell_type','cell_line'):= tstrsplit(elements, ".", fixed=TRUE)
  ][
    ,odds_ratio:=ifelse(p<0.05,odds_ratio,1)
]%>%make_matrix('odds_ratio')%>%log2()
  
  
test_numb_snps = copy(snps_counts_celltypes[c(1:2)])%>%lapply(function(x)x=x[,c(1:4)])
test_numb_snps = Map(mutate,test_numb_snps,ancestry = names(test_numb_snps))%>%rbindlist()%>%make_matrix('numb_snsp_chromstate_celltype')%>%log10()


## color legends for heatmap
# pop
pop_color = my_palette_pop[c(1,3)]
names(pop_color) = c('Denisova','Neanderthal')

enrich_colors=colorRamp2(c(-2,-1,0,1,2), colors=c('blue3','royalblue1','white','red2','red3'))


enrichment_heatmap=function(enrichmatrix,numbsnpsmatrix){
  
  x = Heatmap(enrichmatrix, border = T,
            rect_gp = gpar(type = "none"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.circle(x = x, y = y, r = abs(numbsnpsmatrix[i, j])/25, 
                          gp = gpar(fill = enrich_colors(enrichmatrix[i, j]), col = 'black'))
              
            },
            row_dend_reorder = T,
            row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
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
                              col=pop_color),
              chrom_state=anno_simple(
                gsub("^.*\\.", "", colnames(enrichmatrix)),
                border=T,
                height = unit(2,'cm'),
                col=chromHMM_palette),
              show_annotation_name = F),
            left_annotation = HeatmapAnnotation(which='row',
                                                width = unit(6,'cm'),
                                                Cells = anno_simple(
                                                  gsub("^.*\\.", "",rownames(enrichmatrix)),
                                                  col= tissue_colors),
                                                  show_annotation_name = F))
  return(x)
}

pdf(paste(plot_dir,'enrichment/heatmap_asnp_oddsratio.pdf',sep=''),width=25,height = 50)
enrichment_heatmap(test_enrichment,test_numb_snps)
dev.off()


## write spreadsheet with OR results

table_OR_chromstate = copy(asnps_odds_ratio)[
  ,c(5,2:4,1,8,9,7)
  ]%>%setnames(old='elements',new='chrom_state')

table_OR_chromstate = table_OR_chromstate[order(as.factor(readr::parse_number(table_OR_chromstate$chrom_state))),]

table_OR_chromstate_celltype = rbind(denisovan_odds_ratio_celltype,neanderthal_odds_ratio_celltype)[
  ,c('chrom_state','cell_type','cell_line'):= tstrsplit(elements, ".", fixed=TRUE)
  ][
    ,c(10:12,2:4,1,8,9,7)
]
table_OR_chromstate_celltype = table_OR_chromstate_celltype[order(as.factor(readr::parse_number(table_OR_chromstate_celltype$chrom_state))),]


spreadsheet = list(table_OR_chromstate_celltype,table_OR_chromstate_celltype)
names(spreadsheet) = c('OR chrom state','OR chrom state cell type')

write.xlsx(spreadsheet,paste(table_dir,'Supp_Table_OR_chromstate_pvals.xlsx',sep=''),append=T)

