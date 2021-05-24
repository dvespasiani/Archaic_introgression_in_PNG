library(edgeR)
library(data.table)
library(dplyr)
library(plyr)
library(NineteenEightyR)
library(RColorBrewer)
library(biomaRt)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(viridis)
library(gplots)
library(circlize)

setwd('/Users/dvespasiani/Desktop/Paper_1/')

options(scipen=999)
# # Set output directory and create it if it does not exist:
outputdir="./DE_genes_plots/"

if (file.exists(outputdir) == FALSE){
  dir.create(outputdir)
}
# Load log CPM matrix and y object:
# lcpm
load("./Natri_H_et_al_2020_data/vDup.Rda")
# y DGE list object
load("./Natri_H_et_al_2020_data/voomDupEfit.Rda")

# y DGE list object
load("./Natri_H_et_al_2020_data/indoRNA.read_counts.TMM.filtered.Rda")

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45

y <- calcNormFactors(y, method="TMM")

# assign covariate names
# subtract variables we don't need
subtract=c("group", "norm.factors", "samples")
# get index of unwanted variables
subtract=which(colnames(y$samples) %in% subtract)
covariate.names = colnames(y$samples)[-subtract]
for (name in covariate.names){
  assign(name, y$samples[[paste0(name)]])
}

for (name in c("Age","RIN","lib.size")){
  assign(name, cut(as.numeric(as.character(y$samples[[paste0(name)]])), breaks=5))
}

# assign names to covariate names so you can grab individual elements by name
names(covariate.names)=covariate.names

# assign factor variables
factorVariables=c(colnames(Filter(is.factor,y$samples))[which(colnames(Filter(is.factor,y$samples)) %in% covariate.names)], "Age", "lib.size", "RIN")
numericVariables=colnames(Filter(is.numeric,y$samples))[which(colnames(Filter(is.numeric,y$samples)) %in% covariate.names)] %>% subset(., !(. %in% factorVariables))

# load genes
signif_go_target_genes=as.character(list.files('./target_genes/',recursive = F,full.names = T)) %>% 
  lapply(function(x)fread(x,sep='\t',header=T,select = 6,col.names = 'SYMBOL'))

# nonarch_targets=fread('./top_25_genes_go/non_archaics/png',sep=' ',header=F,col.names='SYMBOL')

# targets=list(arch_targets[[1]],arch_targets[[2]],arch_targets[[3]],nonarch_targets)


# load irene's file and fetch only the DE genes
mtw_kor_de=fread('./Natri_H_et_al_2020_data/topTable.voomNoNorm.tmm.filtered.dup_corrected.MTW-KOR.txt',sep=' ',header = F,drop='V2',
                 col.names = c("genes","logFC","AveExpr", "t", "P.Value","adj.P.Val", "B"))[
                   adj.P.Val<0.05
                   ]


# kor_mtw_smb=rbind(mtw_kor_de,smb_kor_de)

# number of DE genes identified across islands
tot_number_de_genes=copy(mtw_kor_de)[ adj.P.Val<=0.05][,'genes' ] %>% unique() %>% nrow()

# get ensembl ids
mart=useMart('ensembl',dataset="hsapiens_gene_ensembl") 

ensemble_ids=function(x){
  getBM(attributes=c('hgnc_symbol','ensembl_gene_id')
        ,mart = mart,
        filters='hgnc_symbol',
        values=x$SYMBOL)
}
target_ids=lapply(signif_go_target_genes,function(x)ensemble_ids(x))


target_de_ids_kor=lapply(target_ids,function(x)
  x=inner_join(x,mtw_kor_de,by=c('ensembl_gene_id'='genes')) %>% 
    as.data.table()
)

all_de_genes=copy(target_de_ids_kor) %>% 
  lapply(function(x)
    x=x[,c(1:2)] %>% unique())
                  
assign_names=function(x){
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}
all_de_genes=assign_names(all_de_genes)

# # ## write table
# lapply(names(all_de_genes),function(x, all_de_genes)
#   write.table(all_de_genes[[x]], paste(x, "", sep = " "),quote = F,row.names = F),all_de_genes)


arch_target_de_ids_kor=target_de_ids_kor[c(1:2)] %>% 
  lapply(function(x)
  x=x[
    ,c(1,2)
    ] %>% unique()
)
# group the archaics and plot them whereas list as table the naSNPs de genes
# ambig_de_genes=target_de_ids_kor[[1]]
deni_de_genes=target_de_ids_kor[[1]]
nean_de_genes=target_de_ids_kor[[2]]

topGenes=deni_de_genes$hgnc_symbol
topEnsembl=deni_de_genes$ensembl_gene_id

# Let's see how the expression levels of all of the significantly DE genes in population comparisons with Mappi are distributed within each island. First, assign our top genes and ensembl IDs to variables
topGenes=deni_de_genes[hgnc_symbol%like%'OAS']$hgnc_symbol
topEnsembl=deni_de_genes[hgnc_symbol%like%'OAS']$ensembl_gene_id

# To visualise distributions, we'll be making violin plots using ggpubr which needs p-value labels. Let's go ahead and make a matrix to input this into ggpubr
# first set up matrix
topGenes.pvalue=matrix(nrow=length(topEnsembl), ncol=ncol(voomDupEfit))
rownames(topGenes.pvalue)=topEnsembl
colnames(topGenes.pvalue)=colnames(voomDupEfit)
for (i in 1:ncol(voomDupEfit)){
  topTable <- topTable(voomDupEfit, coef=i, n=Inf) # get significant genes over a logFC of 1 for all Island comparisons
  for(j in topEnsembl){
    topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]# input the adjusted p.value for each gene
  }
}


# make pvalues into scientific notation with max 3 digits
topGenes.pvalue=round(topGenes.pvalue,5)
# formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# # convert e notation to base 10 notation
# topGenes.pvalue=sub("e", " x 10^ ", topGenes.pvalue)


# reset ensemble row names to gene symbols
rownames(vDup$E)=vDup$genes$SYMBOL

# We can make the violin plots using ggpubr
pdf(paste0(outputdir,"TopGenes_ggboxplot_Island.pdf"), height=8, width=10)
counter=0
for(ensembl in topEnsembl){
  counter=counter+1
  gene.df <- data.frame(vDup$E[which(vDup$genes$ENSEMBL==ensembl),],Island)
  colnames(gene.df)=c("CPM", "Island")
  annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"),
                              end=c("Mentawai","West Papua","West Papua"), 
                              y=c(max(gene.df[,1]+4),
                                  max(gene.df[,1]+5),max(gene.df[,1]+6)), 
                              label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
  print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island",
                 add=c("jitter","boxplot"), 
                 main=topGenes[counter],
                 palette=c('plum2','olivedrab2','lightskyblue2'),
                 add.params = list(width=0.05),
                 ylab='\n Log2-CPM \n',xlab = ' ') %>% 
          ggpar(legend = 'bottom',ticks = F,font.xtickslab=c(0)) +
          geom_signif(data=annotation_df,
                      aes(xmin=start, xmax=end, annotations=label, y_position=y),
                      textsize = 5, vjust = -0.2,manual=TRUE) + 
          ylim(0, max(gene.df[,1])+7)
  )
}
dev.off()
