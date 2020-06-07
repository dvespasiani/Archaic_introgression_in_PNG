### 
library(data.table);library(magrittr);library(dplyr)
library(tidyr);library(openxlsx)
library(wesanderson);library(RColorBrewer)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

motif_input_dir='./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/'
motif_cluster_dir='./Motifbreak/Motifs_clusters/'
merging_keys=c('seqnames','start','end')

setDTthreads(8)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')


read_tfbs_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
        ,c('start', 'end') :=NULL][
          ,start := snpPos][
            ,end := start +1][
              ,snpPos := NULL]
    )
  x=x[c(2:4)]# remove ambig
  
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}


hocomoco=read_tfbs_snps(paste(motif_input_dir,'hocomoco',sep=''))
jaspar=read_tfbs_snps(paste(motif_input_dir,'jaspar2018',sep=''))
lapply(hocomoco,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
lapply(jaspar,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

# ## motif clusters
motif=fread(paste0(motif_cluster_dir,'motifs_clusters',sep=''),sep=' ',header = T)[
  ,geneSymbol:=ifelse(Motif%like%'H11MO',gsub("\\_HUMAN.*","",Motif),gsub("\\_MA.*","",Motif))
  ][
    ,geneSymbol:=toupper(geneSymbol)
    ][
      ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub('::', '+', geneSymbol)
          ]

tf_cluster=function(x){
  # df=rbind(jaspar_snps_motifbreak,encode_snps_motifbreak,hocomoco_snps_motifbreak)
  # df$motifPos=vapply(df$motifPos, function(x) paste(x, collapse = ","), character(1L))
  df=copy(x)
  df=lapply(df,function(y)y=y[
    ,geneSymbol:=toupper(geneSymbol)
    ][
      ,geneSymbol:=gsub("\\(VAR.2)", "", geneSymbol)
      ][
        ,geneSymbol:=gsub("\\(VAR.3)", "", geneSymbol)
        ][
          ,geneSymbol:=gsub('::', '+', geneSymbol)
          ] %>%inner_join(motif,by='geneSymbol')%>%as.data.table())

}

hocomoco_cluster=tf_cluster(hocomoco)
lapply(hocomoco_cluster,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
lapply(hocomoco,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())

jaspar_cluster=tf_cluster(jaspar)
lapply(jaspar_cluster,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())
lapply(jaspar,function(x)x[,c('seqnames','start','end')] %>% unique() %>% nrow())


combine=function(x,y){
  df=copy(x)
  df1=copy(y)
  df_combined=purrr::map2(df,df1,rbind) %>% lapply(function(a)a=a[,c('seqnames','start','end','Name')] %>% unique())
  return(df_combined)
  
}

tfbs_snps=combine(hocomoco_cluster,jaspar_cluster)

read_active_states=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>%
    lapply(function(y)
      fread(y,sep=' ',header = T,select=c('seqnames','start','end','chrom_state','cell_type','cell_line','all_freq'))[
        chrom_state%in%c('1_TssA','2_TssAFlnk','3_TxFlnk','4_Tx','5_TxWk',"6_EnhG","7_Enh")
        ][
          ,chrom_state:=NULL
          ] %>% unique()
    )
}

states=read_active_states('./Chromatin_states/SNPs_chromHMM_annotated/new_set') # remove the temp directory once the old script is optimized

tfbs_states=purrr::map2(tfbs_snps,states,inner_join,by=merging_keys)
tfbs_states=Map(mutate,tfbs_states,'pop'=names(tfbs_states)) %>% lapply(function(x)setDT(x))
tfbs_states=lapply(tfbs_states,function(x)x=x[
  ,numbsnps_cluster_celltype:=.N,by=.(cell_type,Name)
][,c('Name','cell_type','cell_line','pop','numbsnps_cluster_celltype')] %>% unique())


library(tidyr)
cell_types=copy(tfbs_states[[1]])[,cell_type] %>% unique()
clusters=copy(motif)[,Name] %>% unique()

reference=crossing(cell_types, clusters) %>% as.data.table()
reference=reference[
  ,cluster_cell_type:=paste(clusters,cell_types,sep='.')
][
  ,'cluster_cell_type'
] %>% unique()


test_deni=copy(tfbs_states[[1]])[
  ,cluster_cell_type:=paste(Name,cell_type,sep='.')
][
    ,c('cell_type','cluster_cell_type','numbsnps_cluster_celltype')
    ] %>% unique()


reference_deni=test_deni[reference,on='cluster_cell_type']
reference_deni=reference_deni[
  ,cell_type:=ifelse(is.na(cell_type),gsub(".*\\.","",cluster_cell_type),cell_type)
  ][
    ,numbsnps_cluster_celltype:=ifelse(is.na(numbsnps_cluster_celltype),0,numbsnps_cluster_celltype)
    ][
      ,Cluster:=gsub('\\..*', '', cluster_cell_type)
      ][
        ,cluster_cell_type:=NULL
        ]


combined_references=dcast(reference_deni,cell_type~Cluster,value.var = 'numbsnps_cluster_celltype') 

rownames(combined_references)=combined_references$cell_type
combined_references=combined_references[,-1] %>% as.matrix() 

# combined_references=scale(combined_references)
# combined_references[is.na(combined_references)]=0

# # combined_references[is.na(combined_references)]=1
# 
# 
# 
# 
# test_png=copy(tfbs_states[[2]])[
#   ,cluster_cell_type:=paste(Name,cell_type,sep='.')
#   ][
#     ,pop_celltype:=paste(pop,cell_type,sep='.')
#     ][
#       ,c('pop_celltype','cluster_cell_type','numbsnps_cluster_celltype')
#       ] %>% unique()
# 
# 
# 
# 
# reference_png=test_png[reference,on='cluster_cell_type']
# reference_png=reference_png[
#   ,pop_celltype:=ifelse(is.na(pop_celltype),paste('nean',gsub(".*\\.","",cluster_cell_type),sep='.'),pop_celltype)
#   ][
#     ,numbsnps_cluster_celltype:=ifelse(is.na(numbsnps_cluster_celltype),0,numbsnps_cluster_celltype)
#     ][
#       ,Cluster:=gsub('\\..*', '', cluster_cell_type)
#       ][
#         ,cluster_cell_type:=NULL
#         ]
# 
# 
# combined_references=rbind(reference_deni,reference_png)
# 
# combined_references=dcast(combined_references,pop_celltype~Cluster,value.var = 'numbsnps_cluster_celltype') 
# 
# rownames(combined_references)=combined_references$pop_celltype
# combined_references=combined_references[,-1] %>% as.matrix() %>% na.omit()
# # combined_references[is.na(combined_references)]=1

library(vegan)
test=vegdist(combined_references, method="jaccard")
dist.jac <- as.dist(test)
mds <- as.data.frame(cmdscale(dist.jac))
mds$names <- rownames(mds)

pdf('/home/dvespasiani/test_mds_jaccard.pdf',width = 10,height = 10)
ggplot(mds, aes(V1, V2, label=names)) + 
  # geom_point(aes(colour=factor(cut)), size=2) +
  geom_text(aes(colour=factor(names)), check_overlap = FALSE, size=2.3,
            hjust = "center", vjust = "bottom", nudge_x = 0.005, nudge_y = 0.02) +
  # scale_color_manual(values=c("red", "blue", "green"), guide_legend(title="clusters")) +
  xlab("") + ylab("") +  theme(legend.position = 'none')
dev.off()

# combined_references=scale(combined_references)
# combined_references[is.na(combined_references)]=0

# library(factoextra)
# test=distance(combined_references, method = "jaccard")
# head(test)



# A <- matrix(c(2,5,2,1,0,0,0,0,1,0,0,0,0,1,3,5,6,0,0,1,0,0,0,2,0,0,1,2,7,2,4,6,2,5,1,0,0,1,0,0,0,1,0,0,3,5,4,0,0,1,0,0,1,0,0,2,0,3,5,7,3,1,4,0,1,0,0,0,0,2,0,0,0,1,3,4,6,0,0,1), byrow=T, nrow=8, ncol=10)
# colnames(A) <- letters[1:10]
# rownames(A) <- LETTERS[1:8]
# print(A)
# #weighted jaccard similarity matrix setup
# sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
# rownames(sim.jac) <- rownames(A)
# colnames(sim.jac) <- rownames(A)
# 
# #weighted jaccard function
# # pairs <- t(combn(1:nrow(combined_references), 2))
# for (i in 1:nrow(pairs)){
#   num <- sum(sapply(1:ncol(combined_references), function(x)(min(combined_references[pairs[i,1],x],combined_references[pairs[i,2],x]))))
#   den <- sum(sapply(1:ncol(combined_references), function(x)(max(combined_references[pairs[i,1],x],combined_references[pairs[i,2],x]))))
#   sim.jac[pairs[i,1],pairs[i,2]] <- num/den
#   sim.jac[pairs[i,2],pairs[i,1]] <- num/den  
# }
# sim.jac[which(is.na(sim.jac))] <- 0
# diag(sim.jac) <- 1
# 
# 
# 
# 
# 
# #weighted jaccard similarity matrix setup
# sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
# rownames(sim.jac) <- rownames(A)
# colnames(sim.jac) <- rownames(A)
# 
# #weighted jaccard function
# pairs <- t(combn(1:nrow(A), 2))
# for (i in 1:nrow(pairs)){
#   num <- sum(sapply(1:ncol(A), function(x)(min(A[pairs[i,1],x],A[pairs[i,2],x]))))
#   den <- sum(sapply(1:ncol(A), function(x)(max(A[pairs[i,1],x],A[pairs[i,2],x]))))
#   sim.jac[pairs[i,1],pairs[i,2]] <- num/den
#   sim.jac[pairs[i,2],pairs[i,1]] <- num/den  
# }
# sim.jac[which(is.na(sim.jac))] <- 0
# diag(sim.jac) <- 1
# 
# #weighted jaccard distance
# dist.jac <- 1-sim.jac
# 
# get_dist

test_matrix=combined_references[,-1] %>% as.matrix()
rownames(test_matrix)=combined_references$pop_celltype


pairs <- t(combn(1:nrow(test_matrix), 2))
for (i in 1:nrow(pairs)){
  num <- sum(sapply(1:ncol(test_matrix), function(x)(min(test_matrix[pairs[i,1],x],test_matrix[pairs[i,2],x]))))
  den <- sum(sapply(1:ncol(test_matrix), function(x)(max(test_matrix[pairs[i,1],x],test_matrix[pairs[i,2],x]))))
  sim.jac[pairs[i,1],pairs[i,2]] <- num/den
  sim.jac[pairs[i,2],pairs[i,1]] <- num/den  
}
sim.jac[which(is.na(sim.jac))] <- 0
diag(sim.jac) <- 1




test_deni=copy(tfbs_states[[1]])[
  ,deni:=as.numeric(1)
][
  ,c('Cluster','cell_type','deni')
] %>% unique()

test_deni=test_deni[
  ,cluster_celltype:=paste(Cluster,cell_type,sep='.')
]

test_png=copy(tfbs_states[[3]])[
  ,png:=as.numeric(1)
  ][
    ,c('Cluster','cell_type','png')
    ] %>% unique()




motif_deni=full_join(motif,test_deni,by='Cluster') %>% as.data.table()
test_png=copy(tfbs_states[[3]])
test=rbind(test_deni,test_png)
test=test[
  ,'cell_line':=NULL
]

test=dcast(test,Cluster~cell_type)

# test=test[
#   ,deni:=ifelse(is.na(pop.x),0,1)
# ][
#   ,png:=ifelse(is.na(pop.y),0,1)
#   ][
#     ,c('Cluster','cell_type','cell_line','deni','png')
#   ] %>% unique()


test_adip_deni=copy(test)[cell_line=='Adipose' & pop.x=='denisova'][,c('Cluster','pop.x')] %>% unique()
test_adip_png=copy(test)[cell_line=='Adipose' & pop.y=='png'][,c('Cluster','pop.y')] %>% unique()

motifs_clusters=copy(motif)[,'Cluster'] %>% unique()


check=full_join(motifs_clusters,test_adip_deni,by='Cluster') %>% 
  full_join(test_adip_png,by='Cluster') %>% as.data.table()


test=test[
  ,numb:=.N,by=.(Cluster,cell_type)
][
  ,deni:=ifelse(numb=='1' & pop=='denisova',1,0)
]



