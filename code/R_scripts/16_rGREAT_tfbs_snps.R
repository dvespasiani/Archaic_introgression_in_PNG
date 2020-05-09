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

setDTthreads(8)

target_genes_output_dir='./Motifbreak/GREAT_GO_terms/target_genes/'
de_target_genes_output_dir='./Motifbreak/GREAT_GO_terms/DE_genes_TFBS_regulated/'
top25_go_genes='./Motifbreak/GREAT_GO_terms/top_25_genes_go/'

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

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
  # x=x[c(2:4)]
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') # rm temp dir

tfbs_files=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T,
            drop=c('Refpvalue','Altpvalue','dataSource'))
    )
  x=lapply(x,function(y)
        y=y[
          ,duplicated_effect:=.N,by=.(seqnames,start,end,geneSymbol)
          ]%>% group_by(seqnames,start,end,geneSymbol) %>%
          filter(duplicated_effect==max(duplicated_effect)) %>% as.data.table()
      )
  x=lapply(x,function(y)
          y=y[
            effect%in%'strong'
                    ][
                      ,c('seqnames','start','end','strand','geneSymbol')
                      ] %>% unique()
    )
  x=x[c(2:4)] # remove ambiguous
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  
  x=purrr::map2(x,states,merge,by=c('seqnames','start','end')) # filter snps for those active in t and b cells
  x=lapply(x,function(y)y=as.data.table(y)[all_freq>=0.2][
    ,'all_freq':=NULL
    ])
  return(x)
}

tfbs=tfbs_files('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/combined')

## GREAT enrichments
background=function(x){
  df=copy(x) %>% rbindlist()
    df=df[,geneSymbol:=NULL] %>% unique()%>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
}
deni_background=background(tfbs[c(1,3)])
nean_background=background(tfbs[c(2,3)])
png_background=background(tfbs)

great_enrichment=function(test,background){
  test=copy(test)
  test=test[,geneSymbol:=NULL] %>% unique()
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

pop_names=c('denisova','neandertal','png')
  names(df_list)=pop_names
  for (i in seq_along(df_list)){
    assign(pop_names[i],df_list[[i]],.GlobalEnv)}
  return(df_list)
}

snps_go_enrichment=GO_tables(snps_great)
# 
# write.xlsx(snps_go_enrichment[c(1:2)],'~/pvalue_tables/Supp_Table_9_GO_enriched_terms.xlsx',append=T)

## plot go bp first 25 terms for simplicity of visualization
go_bp=Map(mutate,snps_go_enrichment,'pop'=names(snps_go_enrichment))
go_bp=lapply(go_bp,function(x)
  x=x[1:25,c(1,2,13,14)] %>% as.data.table())

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
      # legend.key=element_blank(),    
      # legend.key.size=unit(1, "cm"),      
      # legend.text=element_text(),
      panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      title=element_text()) +
    coord_flip()
  
}

deni_go_plot=go_enrichment_plot(go_bp[[1]],'#C99E10')
pdf('/home/dvespasiani/go_plot/deni_go_plot.pdf',width=10,height = 7)
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
    Hyper_Adjp_BH<0.01
  ][,1] %>% unique()
)

signif_goterms=lapply(significant_go_terms,function(x)x=x$ID)
# tfbs_target_genes=lapply(snps_gene_targets,function(x)x=unique(x$gene))


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 

genes_assoc_sign_gos=lapply(signif_goterms,function(x)
  x=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'go', values=x, mart = ensembl))
genes_assoc_sign_gos=lapply(genes_assoc_sign_gos,function(x)setDT(x))


target_genes_assoc_sign_gos=purrr::map2(snps_gene_targets,genes_assoc_sign_gos,inner_join,by=c('gene'='hgnc_symbol'))
target_genes_assoc_sign_gos=lapply(target_genes_assoc_sign_gos,function(x)setDT(x))


## get deni genes for cytoscape
deni_genes_for_cytoscape=copy(snps_gene_targets[[1]])
deni_genes_for_cytoscape=getBM(attributes=c('hgnc_symbol','go_id'),filters = 'hgnc_symbol', values=deni_genes_for_cytoscape$gene, mart = ensembl)
deni_top25_go=copy(snps_go_enrichment[[1]][c(1:25),])
deni_genes_for_cytoscape=semi_join(deni_genes_for_cytoscape,deni_top25_go,by=c('go_id'='ID'))
deni_genes_for_cytoscape=deni_genes_for_cytoscape[,1] %>% unique()

write.table(deni_genes_for_cytoscape,paste(top25_go_genes,'denisova_top25_genes',sep=''),col.names = F,quote = F,row.names = F)


## see if any of these target genes is DE between Mentawai and Koorowai
## then look at the fraction of these snps that segregates also in W indonesia
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 

# target_ensembl_ids=lapply(snps_gene_targets,function(x)
# x=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'hgnc_symbol', values=x$gene, mart = ensembl))

# snps_gene_targets=purrr::map2(snps_gene_targets,target_ensembl_ids,inner_join,by=c('gene'='hgnc_symbol'))
# snps_gene_targets=lapply(snps_gene_targets,function(x)setDT(x))


mtw_kor_de=fread('./DE_genes/DE_genes_MTW_KOR.txt',sep=' ',header = F,drop='V2',
                 col.names = c("genes","logFC","AveExpr", "t", "P.Value","adj.P.Val", "B"))[
                   adj.P.Val<0.05
                   ]
tfbs_DE_targets=function(x){
  x=inner_join(x,mtw_kor_de,by=c('ensembl_gene_id'='genes')) %>% 
    as.data.table()
  x=x[,c(1:3,6:8)] %>% unique()
  
}

DE_targets=lapply(target_genes_assoc_sign_gos,function(x)tfbs_DE_targets(x))

## get now CADD scores, alleles in PNG and world wide frequencies
ooa_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      y=fread(y,header = T,sep=' ')
    )
  x=lapply(x,function(y)
    y=y[,c('seqnames','start','end','ref','alt','ANC',
           'POP_ARCH_REF','POP_ARCH_ALT','POP_NOTARCH_REF','POP_NOTARCH_ALT',
           'CADD_RAW','CADD_PHRED','all_frequency','Existing_variation')
        ]%>% unique()
  )
  x=x[c(2:4)]
  return(x)
}

OoA_snps=ooa_snps('./OoA_snps')

DE_targets=purrr:::map2(DE_targets,OoA_snps,inner_join,by=c('seqnames','start','end')) %>% 
  lapply(function(x)setDT(x))


de_target_genes=copy(DE_targets) %>% lapply(function(x)x=x[,c(4,6)] %>% unique())

de_filenames=paste0(de_target_genes_output_dir,names(de_target_genes),sep='')
mapply(write.table,de_target_genes, file = de_filenames,col.names = T, row.names = F, sep = " ", quote = F)


# deni_de_genes_go=de_go_genes(snps_gene_targets[[1]]) # 95
# nean_de_genes_go=de_go_genes(nean_genes_go) #92
# png_de_genes_go=de_go_genes(png_genes_go) #649

## Fetch tfbs-disrupting snps associated with these DE genes
# de_genes=list(ambig_de_genes_go,deni_de_genes_go,nean_de_genes_go,png_de_genes_go)
# 
# 
# target_genes_no_na=lapply(target_genes_no_na,function(x)
#   x=copy(x)[
#     ,c(1:3,6,7)
#     ])
# 
# tfbs_target_de_genes=purrr::map2(target_genes_no_na,de_genes,inner_join,by=c('gene'='hgnc_symbol')) %>% 
#   purrr::map2(tfbs,inner_join,by=c('seqnames','start','end')) %>% 
#   lapply(function(x)setDT(x))


## see how many are within W indonesia
denisova_in_w_indo=fread('../Original_files_per_region/Denisova/deni_w_indo',sep='\t',header = T,drop=c(11:14))
neandertal_in_w_indo=fread('../Original_files_per_region/Neandertal/nean_w_indo',sep='\t',header = T,drop=c(11:14))

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


deni_shared=shared_snps(DE_targets[[1]],denisova_in_w_indo) # 78/284 shared  between w indo and png
nean_shared=shared_snps(DE_targets[[2]],neandertal_in_w_indo) # 158/224 shared between w indo and png
png_shared=shared_snps(DE_targets[[3]],windo_nasnps) # 252/512 shared between w indo and png


# barplot alleles shared between png and w indo
windo_png=function(x,y){
  x=copy(x)%>% 
    anti_join(y[,c(1:3)],by=c('seqnames','start','end')) %>%
    as.data.table()
  x=x[
    ,'distribution':='png'
    ][
      ,c('seqnames','start','end','distribution')
      ] %>% unique()
  x=rbind(x,y)
  x=x[
    ,'tot_snps':=.N
    ][
      ,'tot_snps_per_condition':=.N,by=.(distribution)
      ][
        ,'distribution_fraction':=tot_snps_per_condition/tot_snps
        ][
          ,c('distribution','distribution_fraction')
          ] %>% unique()
}


windo_png_neandertal=windo_png(DE_targets[[2]],nean_shared)[,'pop':='Neandertal']
windo_png_denisova=windo_png(DE_targets[[1]],deni_shared)[,'pop':='Denisova']

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



pdf('/home/dvespasiani/tfbs_plots/TFBS_DE_genes_ISEA_distribution.pdf',width = 7,height = 7)
shared_plot
dev.off()

# ## cadd scores of these TFBS-snps targeting DE genes
# cadd_scores_plot=function(x){
#   y=copy(x)
#   y=Map(mutate,y,'pop'=names(y)) %>% rbindlist()
#   y=y[,c(1:3,14,18)] %>% unique()
#   
#   my_palette_pop=c('#C99E10', # denisova
#                    '#9B4F0F', #neandertal
#                    '#1E656D' #png
#   )
#   
#   names(my_palette_pop)= levels(as.factor(y$pop))
#   colScale=scale_fill_manual(name= " ",values = my_palette_pop,labels = c('Denisova','Neandertal','PNG'))
#   
#   plot=ggplot(y,aes(x=pop,y=CADD_RAW,fill=pop))+
#     colScale+
#     xlab(' ')+ylab('\n CADD raw score \n ')+
#     geom_violin(trim=F,scale = "width")+
#     geom_boxplot(width=.1, position =  position_dodge(width = 0.4),notch = T)+
#     stat_compare_means(method = "wilcox.test",label = "p.signif",
#                        ref.group = "png",hide.ns =T,
#                        label.y =max(y$CADD_RAW)+0.00003,
#                        size=12)+
#     theme(axis.text.x = element_blank(),
#           axis.text.y = element_text(size=30),
#           axis.title.y = element_text(size=40,hjust=0.5),
#           axis.title.x = element_blank(),
#           axis.ticks.x = element_blank(),
#          plot.title = element_text(hjust = 0.5,size=15),
#           plot.subtitle = element_text(hjust= 0.5),
#           panel.background =element_rect(fill = 'white', colour = 'black'),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           legend.text = element_text(size = 28),
#           legend.title = element_text(size = 28),
#           legend.margin = margin(c(0.5, 2, 8, 25)),
#           legend.spacing.x = unit(0.5, 'cm'),
#           axis.line = element_line(color = "black",
#                                    size = 0.7, linetype = "solid"))
#   return(plot)
# }
# 
# pdf('/home/dvespasiani/tfbs_plots/TFBS_de_genes_cadd_plot.pdf',width = 20,height = 20)
# cadd_scores_plot(DE_targets)
# dev.off()

## get OAS genes TFBS-snps
OAS_snps=copy(target_genes_assoc_sign_gos[[1]])[gene%in%c('OAS2','OAS3')] %>% inner_join(tfbs[[1]],by=c('seqnames','start','end'))
OAS_snps_mtw_kor=inner_join(OAS_snps,denisova_in_w_indo,by=c('seqnames'='CHR','start'='FROM','end'='TO')) %>% as.data.table()

## see which TFBS they disrupt

# 
# shared_snps=function(x,archaic){
#   x=x[,c('seqnames','start','end','geneSymbol','gene')]  %>% 
#     inner_join(archaic,by=c('seqnames'='CHR','start'='FROM','end'='TO')) %>% 
#     filter(`POP_ARCH_REF`+`POP_ARCH_ALT`>1) %>% # removes singletons
#     as.data.table()
#   x=x[
#       ,c(1:3)
#       ][
#         ,'distribution':='shared'
#         ] %>% unique()
# }
# 
# deni_shared=shared_snps(tfbs_target_de_genes[[2]],denisova_in_w_indo) # 84/348  between w indo and png
# nean_shared=shared_snps(tfbs_target_de_genes[[3]],neandertal_in_w_indo) # 252/322 shared between w indo and png
# 
# ambig_shared_nean=shared_snps(tfbs_target_de_genes[[1]],neandertal_in_w_indo)
# ambig_shared_deni=shared_snps(tfbs_target_de_genes[[1]],denisova_in_w_indo)           
# 
# ambig_shared=rbind(ambig_shared_nean,ambig_shared_deni) %>% unique() # 68/123
# 
# # barplot 
# windo_png=function(x,y){
#   x=copy(x)%>% 
#   anti_join(y[,c(1:3)],by=c('seqnames','start','end')) %>%
#     as.data.table()
#   x=x[
#       ,'distribution':='png'
#     ][
#   ,c('seqnames','start','end','distribution')
# ] %>% unique()
# x=rbind(x,y)
# x=x[
#   ,'tot_snps':=.N
# ][
#   ,'tot_snps_per_condition':=.N,by=.(distribution)
# ][
#   ,'distribution_fraction':=tot_snps_per_condition/tot_snps
# ][
#   ,c('distribution','distribution_fraction')
# ] %>% unique()
# }
# 
#    
# windo_png_neandertal=windo_png(tfbs_target_de_genes[[3]],nean_shared)[
#   ,'pop':='Neandertal'
# ]
# windo_png_denisova=windo_png(tfbs_target_de_genes[[2]],deni_shared)[
#   ,'pop':='Denisova'
# ]
# 
# windo_png_ambiguous=windo_png(tfbs_target_de_genes[[1]],ambig_shared)[
#   ,'pop':='Ambiguous'
#   ]
# 
# windo_png_archaics=rbind(windo_png_denisova,windo_png_neandertal,windo_png_ambiguous)
# 
# shared_plot=ggplot(windo_png_archaics,aes(x=pop,y=distribution_fraction,fill=distribution))+
#   geom_bar(stat='identity')+
#   scale_fill_manual(values=c('lightskyblue2','plum2'),
#                     name=' ',labels=c('Papuan specific','Shared with west Indonesians'))+
#   ylab('\n Fraction of SNPs \n')+ xlab('\n \n')+
# theme(axis.text.x = element_text(hjust = 1,size=28,angle = 40),
#       axis.text.y = element_text(size=28),
#       axis.title.y = element_text(size=30,hjust=0.5),
#       axis.text=element_text(size=28),
#       legend.text = element_text(size = 28),
#       legend.title = element_text(size = 28),
#       legend.margin = margin(c(0.5, 2, 8, 25)),
#       legend.spacing.x = unit(0.5, 'cm'),
#       axis.line = element_line(color = "black", 
#                                size = 0.7, linetype = "solid"))
# 
#   
# 
# pdf('/home/dvespasiani/TFBS_DE_genes_ISEA_distribution.pdf',width = 20,height = 20)
# shared_plot
# dev.off()
# 
# ## papuan specific TFBS-disrupting aSNPs
# papuans_deni=anti_join(tfbs_target_de_genes[[2]],deni_shared,by=c('seqnames','start','end')) %>%
#   as.data.table()
# papuans_oas=copy(papuans_deni)[
#   ,c('seqnames','start','end','REF','ALT','geneSymbol','gene')
# ][gene%in%c('OAS1','OAS2','OAS3')] %>% unique() # 12 snps
# 
# # 3 snps are not present in w indonesians and for the others these people are fixed for the ref allele
# oas_snps=inner_join(papuans_oas,denisova_in_w_indo,by=c('seqnames'='CHR','start'='FROM','end'='TO','REF','ALT')) %>%
#   inner_join(ooa_snps[[2]],by=c('seqnames','start','end','ANC','REF'='ref','ALT'='alt')) %>% 
#   as.data.table() 
# 
# ## possible cool snps are those that disrupt immune specific tfs (e.g. IRF4/TCF12)
# denisova_oas_cool_snps=copy(oas_snps)[
#   geneSymbol%in%'IRF4'
#   ] # rs372433785_chr12_113345315_G_C 
# 
# 
