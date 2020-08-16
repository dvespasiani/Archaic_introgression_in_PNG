# rGREAT tfbs snps strong effect with freq >0.2 in T/B cells
library(dplyr); library(data.table)
library(magrittr); 
library(purrr)
library(rGREAT)
library(GenomicRanges)
library(biomaRt)
library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

target_genes_output_dir='./Motifbreak/GREAT_GO_terms/target_genes/'
de_target_genes_output_dir='./Motifbreak/GREAT_GO_terms/DE_genes_TFBS_regulated/'
top25_go_genes='./Motifbreak/GREAT_GO_terms/top_25_genes_go/'

table_dir='./Results/Tables/'
plot_dir='./Results/Plots/GREAT/'

allele_frequency=0.3


setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')



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


read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq'))[
        chrom_state%in%c('2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
        ][
          cell_line%in%c('TCells','BCells')
          ][
            ,c('chrom_state','cell_line','cell_type'):=NULL
            ] %>% unique()
    )
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') 

read_tfbs=function(x,pwm_effect){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        !dataSource%like%'HOCOMOCOv11-secondary-D'][
          ,end:=start+1
        ][
          ,c('seqnames','start','end','effect')
        ][effect%in%pwm_effect] %>% unique())
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  
  x=purrr::map2(x,states,merge,by=c('seqnames','start','end')) 
  x=lapply(x,function(y)y=as.data.table(y)[,all_freq:=round(all_freq,1)][all_freq>=allele_frequency][
    ,c('all_freq','effect'):=NULL
    ] %>% unique())
  return(x)
}

tfbs=read_tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/new_set/combined','strong')
lapply(tfbs,function(y)y=y[,c(1:3)] %>% unique() %>% nrow())

## GREAT enrichments
background=function(x){
  df=copy(x) %>% rbindlist()%>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
}
deni_background=background(tfbs[c(1,3)])
nean_background=background(tfbs[c(2,3)])
png_background=background(tfbs)

great_enrichment=function(test,background){
  test=copy(test) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  result=submitGreatJob(gr=test,
                   bg=background,
                   species = "hg19",
                   rule= "twoClosest", 
                   adv_twoDistance = 1000,# distances are in kb (taken from GTEx observed TFBS cis-eQTLs)
                   includeCuratedRegDoms = T,
                   request_interval=10)

}
deni_great=great_enrichment(tfbs[[1]],deni_background)
nean_great=great_enrichment(tfbs[[2]],nean_background)
png_great=great_enrichment(tfbs[[3]],png_background)

snps_great=list(deni_great,nean_great,png_great)

GO_tables=function(x){
  df=copy(x)
  df=lapply(df,function(y)
    y=getEnrichmentTables(y,ontology='GO Biological Process'))
  
df_list=list(df[[1]][[1]],df[[2]][[1]],df[[3]][[1]])
df_list=lapply(df_list,function(x)x=as.data.table(x)[Hyper_Adjp_BH<=0.01])
pop_names=c('denisova','neandertal','png')
  names(df_list)=pop_names
  for (i in seq_along(df_list)){
    assign(pop_names[i],df_list[[i]],.GlobalEnv)}
  return(df_list)
}

snps_go_enrichment=GO_tables(snps_great)
head(snps_go_enrichment[[1]],30)
head(snps_go_enrichment[[2]],30)

write.xlsx(snps_go_enrichment,paste(table_dir,'Supp_Table_10_GO_enriched_terms.xlsx',sep=''),append=T)

## plot go bp first 25 terms for simplicity of visualization
go_bp=Map(mutate,snps_go_enrichment,'pop'=names(snps_go_enrichment))
go_bp=lapply(go_bp,function(x) x=x[1:30,c(1,2,13,14)] %>% as.data.table())

go_enrichment_plot=function(x,y){
  ggplot(x, aes(x=reorder(name,-log10(Hyper_Adjp_BH)), y=-log10(Hyper_Adjp_BH),fill=pop)) +
    geom_bar(stat = 'identity',position = 'dodge')+
    xlab(" ") +
    ylab("\n -Log10 (P) \n ") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(x$Hyper_Adjp_BH)), by = 2), 1)) +
    scale_fill_manual(name= " ",values=y)+
   theme(
      legend.position='none',
      legend.background=element_rect(),
      axis.text.x=element_text(angle=0, hjust=1.10),
      axis.text.y=element_text(angle=0, vjust=0.8),
      axis.title=element_text(),
      axis.line = element_line(color = "black",size = 0, linetype = "solid"),
     panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      title=element_text()) +
    coord_flip()
  
}

deni_go_plot=go_enrichment_plot(go_bp[[1]],'#C99E10')
pdf(paste(plot_dir,'deni_go_plot.pdf',sep=''),width=10,height = 7)
deni_go_plot
dev.off()

# get list of target genes 
target_genes= function(snps){
  genes=lapply(snps,function(x)x=plotRegionGeneAssociationGraphs(x,type=1,request_interval = 10))
  genes=lapply(genes,function(x)as.data.table(x) %>% na.omit() %>% unique())
  
  pop_names=c('denisova','neandertal','png')
  names(genes)=pop_names
  for (i in seq_along(genes)){
    assign(pop_names[i],genes[[i]],.GlobalEnv)}
  return(genes)
  }

snps_gene_targets=target_genes(snps_great)

filenames=paste0(target_genes_output_dir,names(snps_gene_targets),sep='')
mapply(write.table,snps_gene_targets, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)

## Fetch TFBS-target genes related to the significant GO terms
significant_go_terms=copy(snps_go_enrichment) %>%
  lapply(function(x)
  x=as.data.table(x)[
    Hyper_Adjp_BH<=0.01
  ][,1] %>% unique()
)

signif_goterms=lapply(significant_go_terms,function(x)x=x$ID)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 

genes_assoc_sign_gos=lapply(signif_goterms,function(x)
  x=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'go', values=x, mart = ensembl))
genes_assoc_sign_gos=lapply(genes_assoc_sign_gos,function(x)setDT(x))

target_genes_assoc_sign_gos=purrr::map2(snps_gene_targets,genes_assoc_sign_gos,inner_join,by=c('gene'='hgnc_symbol'))
target_genes_assoc_sign_gos=lapply(target_genes_assoc_sign_gos,function(x)setDT(x))

## get deni genes for cytoscape
all_deni_targets=copy(target_genes_assoc_sign_gos[[1]])
all_deni_targets=all_deni_targets$gene%>%unique()

deni_genes_for_cytoscape=copy(snps_gene_targets[[1]])
deni_genes_for_cytoscape=getBM(attributes=c('hgnc_symbol','go_id'),filters = 'hgnc_symbol', values=deni_genes_for_cytoscape$gene, mart = ensembl)
deni_top25_go=copy(snps_go_enrichment[[1]][c(1:30),])

deni_genes_for_cytoscape=semi_join(deni_genes_for_cytoscape,deni_top25_go,by=c('go_id'='ID'))
deni_genes_for_cytoscape=deni_genes_for_cytoscape[,1] %>% unique()

write.table(deni_genes_for_cytoscape,paste(top25_go_genes,'denisova_top25_genes',sep=''),col.names = F,quote = F,row.names = F)
write.table(all_deni_targets,paste(top25_go_genes,'denisova_all_targets',sep=''),col.names = F,quote = F,row.names = F)
###-------------------------------------------------------------------------------
##     see if any of these target genes is DE between Mentawai and Koorowai
##   then look at the fraction of these snps that segregates also in W indonesia
###-------------------------------------------------------------------------------
mtw_kor_de=fread('./DE_genes/DE_genes_MTW_KOR.txt',sep=' ',header = F,drop='V2',
                 col.names = c("genes","logFC","AveExpr", "t", "P.Value","adj.P.Val", "B"))[
                   adj.P.Val<=0.01
                   ]
tfbs_DE_targets=function(x){
  x=inner_join(x,mtw_kor_de,by=c('ensembl_gene_id'='genes')) %>% 
    as.data.table()
  x=x[,c(1:3,6:8)] %>% unique()
  
}

DE_targets=lapply(target_genes_assoc_sign_gos,function(x)tfbs_DE_targets(x))

## test for enrichment in numb de genes targeted by deni/nean/png vs all human genes
## get all human genes via annotatr
annots = c('hg19_basicgenes')
human_genes = annotatr::build_annotations(genome = 'hg19', annotations = annots) %>%
  as.data.table()
human_genes=human_genes[,'symbol'] %>% na.omit() %>% unique() %>% nrow()

enrichment_pval=function(x,y){
  
  numb_de_genes_per_pop=copy(x)[,'ensembl_gene_id'] %>% unique() %>% nrow() 
  numb_genes_associated_go=copy(y)[,'ensembl_gene_id'] %>% unique() %>% nrow()
  tot_numb_de_genes=copy(mtw_kor_de)[,'genes'] %>% unique() %>% nrow() 
  
  test=phyper(numb_de_genes_per_pop-1,tot_numb_de_genes,human_genes-tot_numb_de_genes,numb_genes_associated_go,lower.tail = F)
  return(test)
}

enrichment_pval(DE_targets[[1]],target_genes_assoc_sign_gos[[1]]) ## 0.0020
enrichment_pval(DE_targets[[2]],target_genes_assoc_sign_gos[[2]]) ## 0.049
enrichment_pval(DE_targets[[3]],target_genes_assoc_sign_gos[[3]]) ## 1.396314e-15

lapply(DE_targets,function(x)x=x[,'gene']%>%unique()%>%nrow())
lapply(target_genes_assoc_sign_gos,function(x)x=x[,'gene']%>%unique()%>%nrow())

## get now alleles in PNG 
DE_targets=purrr:::map2(DE_targets,original_snps,inner_join,by=c('seqnames','start','end')) %>% 
  lapply(function(x)setDT(x))
lapply(DE_targets,function(x)x=x[,'gene'] %>% unique() %>% nrow())

de_target_genes=copy(DE_targets) %>% lapply(function(x)x=x[,c(4,6)] %>% unique())

# de_filenames=paste0(de_target_genes_output_dir,names(de_target_genes),sep='')
# mapply(write.table,de_target_genes, file = de_filenames,col.names = T, row.names = F, sep = " ", quote = F)

## see how many are within W indonesia
denisova_in_w_indo=fread('./Original_files_per_region/Denisova/deni_w_indo.gz',sep='\t',header = T,drop=c(11:14))
neandertal_in_w_indo=fread('./Original_files_per_region/Neandertal/nean_w_indo.gz',sep='\t',header = T,drop=c(11:14))

shared_snps=function(de,windo_snps){
  de=de[,c('seqnames','start','end','gene')]  %>% 
    inner_join(windo_snps,by=c('seqnames'='CHR','start'='FROM','end'='TO')) %>% 
    filter(`POP_ARCH_REF`+`POP_ARCH_ALT`>1) %>% # removes singletons
    as.data.table()
  de=de[
    ,c(1:3)
    ][
      ,'distribution':='shared'
      ] %>% unique()
}


deni_shared=shared_snps(DE_targets[[1]],denisova_in_w_indo) 
nean_shared=shared_snps(DE_targets[[2]],neandertal_in_w_indo)

# barplot alleles shared between png and w indo
windo_png=function(x,y){
  png_specific=copy(x)%>% 
    anti_join(y[,c(1:3)],by=c('seqnames','start','end')) %>%
    as.data.table()
  png_specific=png_specific[
    ,'distribution':='png'
    ][
      ,c('seqnames','start','end','distribution')
      ] %>% unique()
  isea_distribution=rbind(png_specific,y)
  isea_distribution=isea_distribution[
    ,'tot_snps':=.N
    ][
      ,'tot_snps_per_condition':=.N,by=.(distribution)
      ][
        ,'distribution_fraction':=tot_snps_per_condition/tot_snps
        ][
          ,c('distribution','distribution_fraction','tot_snps','tot_snps_per_condition')
          ] %>% unique()
  return(isea_distribution)
}


windo_png_denisova=windo_png(DE_targets[[1]],deni_shared)[,'pop':='Denisova']
windo_png_neandertal=windo_png(DE_targets[[2]],nean_shared)[,'pop':='Neandertal']

windo_png_archaics=rbind(windo_png_neandertal,windo_png_denisova)

shared_plot=ggplot(windo_png_archaics,aes(x=pop,y=distribution_fraction,fill=distribution))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c('lightskyblue2','plum2'),
                    name=' ',labels=c('Papuan specific','Shared'))+
  ylab('\n Fraction of TFBS-SNPs \n')+ xlab('\n \n')+
  theme(axis.text.x = element_text(hjust = 1,angle = 40),
        axis.text.y = element_text(),
        axis.title.y = element_text(hjust=0.5),
        axis.text=element_text(),
        legend.text = element_text(),
        legend.title = element_text(),
        legend.spacing.x = unit(0.5, 'cm'),
        panel.background =element_rect(fill = 'white', size = 1,colour = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0, linetype = "solid"))



pdf(paste(plot_dir,'TFBS_DE_genes_ISEA_distribution.pdf',sep=''),width = 7,height = 7)
shared_plot
dev.off()


## get OAS genes TFBS-snps
OAS_snps=copy(target_genes_assoc_sign_gos[[1]])[gene%in%c('OAS2','OAS3')] %>%  ## only for these genes as they are DE 
  inner_join(tfbs[[1]],by=c('seqnames','start','end'))
OAS_snps_mtw=inner_join(OAS_snps,denisova_in_w_indo,by=c('seqnames'='CHR','start'='FROM','end'='TO')) %>% as.data.table()

OAS_snps_mtw_png=OAS_snps_mtw[original_snps[[1]],on=c('seqnames','start','end','REF',"ALT",'ANC'),nomatch=0] ## freq difference between mtw and krw





