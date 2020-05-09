## tfbs that are also gtex eqtls
library(data.table);library(magrittr);library(dplyr)
library(purrr)
library(GenomicRanges);library(biomaRt);
library(openxlsx)
library(ggplot2);library(wesanderson);library(ggthemes);library(ggfortify)

setDTthreads(8)

setwd('/data/projects/punim0586/dvespasiani/Files/PNG')

# gtex v8 tissue specific cis-eqtls 
cis_eqtls=as.character(list.files('./GTEx/GTEx_Analysis_v8_eQTL/',full.names = T,recursive = F)) %>% 
  lapply(function(x) x=fread(x,sep = '\t',header = T,
          select = c('chr','variant_pos','ref','alt','tss_distance', 'gene_id',
                      'gene_name','gene_chr','gene_start','gene_end','strand',
                      'num_alt_per_site','rs_id_dbSNP151_GRCh38p7',
                      'maf','qval','pval_nominal_threshold','log2_aFC','log2_aFC_lower','log2_aFC_upper'
                      ))
         [
                        qval<=0.05
                      ]
)

assign_names=function(x,list_df){
  gtex_tissues=as.character(list.files(x,full.names = F,recursive = F)) 
  gtex_tissues=gsub("\\..*","",gtex_tissues)
  names(list_df)=gtex_tissues
  for (i in seq_along(x)){
    assign(gtex_tissues[i],x[[i]],.GlobalEnv)}
  return(list_df)
}

cis_eqtls=assign_names('./GTEx/GTEx_Analysis_v8_eQTL/',cis_eqtls)
cis_eqtls=Map(mutate,cis_eqtls,'tissue'=names(cis_eqtls))

# separate liftover for cis-regulatory variant and gene
cis_eqtls_snps=copy(cis_eqtls)

cis_eqtls_snps=lapply(cis_eqtls_snps,function(x)
  x=x[ ,c('chr','variant_pos','ref','alt','rs_id_dbSNP151_GRCh38p7','tissue','gene_id','gene_name')
] %>% makeGRangesFromDataFrame(seqnames.field = 'chr',start.field = 'variant_pos',
                               end.field = 'variant_pos', keep.extra.columns = T))

regulated_genes=copy(cis_eqtls)
regulated_genes=lapply(regulated_genes,function(x)
  x=x[
  ,c('gene_chr','gene_start','gene_end','strand','gene_id','gene_name','maf','qval','pval_nominal_threshold','log2_aFC', 'log2_aFC_lower','log2_aFC_upper','tissue')
]%>% makeGRangesFromDataFrame(seqnames.field = 'gene_chr',start.field = 'gene_start',
                              end.field = 'gene_end',strand.field = 'strand', keep.extra.columns = T)
)


# liftover to change from hg38 to hg19
library(rtracklayer)
liftover_chain=import.chain('./LiftOver/hg38ToHg19.over.chain')

liftover=function(x,y){
  x=liftOver(x,y) 
  x=unlist(x) %>% as.data.table()
}

cis_eqtls_snps_hg19 = lapply(cis_eqtls_snps,function(x)liftover(x,liftover_chain))
cis_eqtls_snps_hg19=lapply(cis_eqtls_snps_hg19,function(x)x=x[,c('width','strand'):=NULL])
regulated_genes_hg19 = lapply(regulated_genes,function(x)liftover(x,liftover_chain))

cis_eqtls_hg19=purrr::map2(cis_eqtls_snps_hg19,regulated_genes_hg19,merge,by=c('gene_id','gene_name','tissue'))
  
cis_eqtls_hg19=lapply(cis_eqtls_hg19,function(x)
  x=x[
   ,.(snp_seqnames=seqnames.x,snp_start=start.x,ref,alt,rs_id_dbSNP151=rs_id_dbSNP151_GRCh38p7,
      gene_seqnames=seqnames.y,gene_start=start.y,gene_end= end.y, width,strand,maf, 
      gene_id,gene_name,tissue,
      qval,pval_nominal_threshold,log2_aFC,log2_aFC_lower,log2_aFC_upper)
  ]
  )
cis_eqtls_hg19=rbindlist(cis_eqtls_hg19)

# read tfbs and merge dfs
tfbs=function(x){
 x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
fread(y,sep=' ',header = T,
            drop=c('Refpvalue','Altpvalue','dataSource')))
 
 x=x[c(2:4)] #remove ambiguous
pop_names=c('denisova','neandertal','png')
names(x)=pop_names
for (i in seq_along(x)){
  assign(pop_names[i],x[[i]],.GlobalEnv)}
return(x)

}
 
snps_tfbs=tfbs('./Motifbreak/Tx_and_CREs/TFBSs_disrupted_10neg5/combined')
lapply(snps_tfbs,function(x)x=x[,c('seqnames','start','end')] %>% unique() %>% nrow())
## ambiguous 5922
## denisova 42689
## neandertal 24292
## png 259386

## SNPS that are eQTLs
tfbs_eqtls=function(x){
  x=lapply(x,function(y)inner_join(y,cis_eqtls_hg19,by=c('seqnames'='snp_seqnames','start'='snp_start', 'REF'='ref')) %>% 
             as.data.table())
  x=Map(mutate,x,'pop'=names(x))
  y=lapply(x,function(y)setDT(y))
}
snps_tfbs_eqtls=tfbs_eqtls(snps_tfbs)

## table tfbs snps that are cis-eQTLs
table_tfbs_eqtls=lapply(snps_tfbs_eqtls,function(x)
  x=x[
    ,c('seqnames','start','end','rs_id_dbSNP151','alt','gene_name','tissue','maf','qval')
  ] %>% setnames(c('seqnames','start','end','rsid','risk_allele','gene_name','tissue','maf','qval'))%>% unique()
    )
# lapply(table_tfbs_eqtls,function(x)x=x[,c(1:3)] %>% unique() %>% nrow())
# xxx ambiguous
# 34 denisova (0.079%)
# 76 neandertal (0.31%)
# 402 png (0.15%)

## write tfbs eqtls
# tfbs_output_dir='./Motifbreak/TFBS_eQTLs/'
# filenames=paste0(tfbs_output_dir,names(table_tfbs_eqtls),sep='')
# mapply(write.table,snps_tfbs_eqtls, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)

# frequency differences between gtex maf and snp freq in png
original_files=function(x){
  x=as.character(list.files(x,full.names = T,recursive = F)) %>% 
    lapply(function(y)
      fread(y,sep=' ',header = T)[
       ,freq_range:=ifelse(all_frequency<=0.05,'low','high') # adjusted  
      ]
    )
  x=x[c(2:4)]
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  x=Map(mutate,x,'pop'=names(x)) 
  x=lapply(x,function(y)y %>% as.data.table())
  return(x)
  
}
snps=original_files('./Grouped_filtered_snps/')

gtex_png_comparison=function(png,gtex){
  png=copy(png)
  gtex=copy(gtex)
  png=lapply(png,function(x)
    x=x[,frequency_type:='png_snpfreq'
        ][
          , snp_freq:=all_frequency
          ][
            ,c('CHR','FROM','TO','snp_freq','pop','frequency_type')
            ] %>% unique())
  
  gtex=lapply(gtex,function(x)
    x=x[,frequency_type :='maf'
        ][
          ,snp_freq:=maf
          ][
            ,c('seqnames','start','end','snp_freq','pop','frequency_type') 
            ] %>% unique()%>% group_by(`seqnames`,`start`,`end`) %>% filter(`snp_freq`==max(snp_freq))%>% 
      as.data.table())
  
  x=purrr::map2(png,gtex,semi_join,by=c('CHR'='seqnames','FROM'='start','TO'='end'))
  y=purrr::map2(gtex,png,semi_join,by=c('seqnames'='CHR','start'='FROM','end'='TO')) 
  x=rbindlist(x)[
    ,c(4:6)
    ]
  y=rbindlist(y)[
    ,c(4:6)
    ]
  z=rbind(x,y)
  return(z)
 
}

snps_eqtl_freq=gtex_png_comparison(snps,snps_tfbs_eqtls)

pvalues=copy(snps_eqtl_freq)
pvalues=pvalues%>% split(as.factor(snps_eqtl_freq$pop)) %>% 
  lapply(function(x)x=x[
    ,p:=wilcox.test(snp_freq ~ frequency_type, data=x,paired=T)$p.val
    ][
      ,statistic:=wilcox.test(snp_freq ~ frequency_type, data=x,paired=T)$statistic
    ][
      ,p.adj:=p.adjust(p,'bonferroni')
    ][
      ,p.signif:= ifelse(`p.adj`<=0.0001,'****',
                         ifelse(`p.adj`>0.0001 &`p.adj`<=0.001,'***',
                                ifelse(`p.adj`>0.001 & `p.adj`<=0.01,'**',
                                       ifelse(`p.adj`>0.01 & `p.adj`<=0.05,'*',' '))))
      ][
        ,c('pop','statistic','p','p.adj','p.signif')
      ] %>% unique()) %>% rbindlist()

## p values
## deni (W=335, p=5.317529e-01,padj= 1.000000e+00)
## nean (W=556, p=2.688617e-06, padj=4.086698e-04)
## png (W=20566, p=1.211664e-17,padj= 9.741775e-15)

eQTL_freq_plot=function(df){
  
  pval=copy(pvalues)[,frequency_type:='png_snpfreq']
  
  
  my_palette=c('#999999', # GTEx V8
               'lightgoldenrod3' #Jacobs et al.2019
  )
  
  names(my_palette)= levels(as.factor(df$frequency_type))
  colScale=scale_fill_manual(name= " ",values = my_palette,labels=c('GTEx V8','Jacobs et al.2019')) 
  
  ggplot(df,aes(x=pop,y=snp_freq,fill=frequency_type))+
    geom_violin(trim=F,scale = "width")+
    geom_boxplot(width=0.2, position =  position_dodge(width = 0.9),outlier.size = 0.3)+
    # geom_jitter(position=position_jitter(0.1))+
    geom_text(data=pval, aes(x=pop, y=1.0001, label=p.signif), col='black', size=7)+
    scale_x_discrete(limit = c("denisova", "neandertal",'png'),labels = c("Denisova","Neandertal",'Papuans'))+
    colScale+
    xlab(' ')+
    ylab('\n eQTL SNPs frequency \n')+
    theme(panel.background =element_rect(fill = 'white', colour = 'white'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(),
          legend.title = element_text(),
          legend.margin = margin(c(0.5, 2, 8, 25)),
          legend.spacing.x = unit(0.5, 'cm'),
          legend.position = 'bottom',
          axis.text.x = element_text(angle = 60,vjust = 0.4,hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(),
          axis.title.y = element_text(hjust=0.5),
          axis.text=element_text(),
          axis.line = element_line(color = "black",size = 0.5, linetype = "solid"))
  
}

pdf('/home/dvespasiani/tfbs_plots/tfbs_eqtls_freq.pdf',width = 7,height =7)
eQTL_freq_plot(snps_eqtl_freq)
dev.off()

## check if tfbs snp is on the same strand of regulated gene and how far they are ~
## median distance for denisova = 1kb for nean 26kb for ambiguous 12kb for png 14kb
tss_distance=copy(snps_tfbs_eqtls)
tss_distance=lapply(tss_distance,function(x)
  x=x[
    ,tss_distance:=gene_start-start
  ][
    ,same_strand:=ifelse(strand.x==strand.y,'yes','no')
  ][,c('seqnames','start','end','tss_distance','pop')] %>% 
    unique())

tss_distance=rbindlist(tss_distance)[
  ,c('tss_distance','pop')
][
  ,tss_distance:=round(log10(tss_distance),1)
] %>% na.omit()


tss_distance_plot=ggplot(tss_distance,aes(x=pop,y=tss_distance,fill=pop))+
  geom_violin(trim=F,scale = "width")+
  geom_boxplot(width=0.2, position =  position_dodge(width = 0.9),outlier.size = 0.3)+
  scale_fill_manual(values=c('#C99E10','#9B4F0F','#1E656D'),name='Population',labels=c('Denisova','Neandertal','PNG'))+
  scale_x_discrete(limit = c("denisova", "neandertal",'png'),labels = c("Denisova","Neandertal",'Papuans'))+
  xlab('\n  \n')+ylab('\n log10 base distance eQTL SNPs - eGene TSS  \n')+
  theme(panel.background =element_rect(fill = 'white', colour = 'white'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(),
        legend.title = element_text(),
        legend.margin = margin(c(0.5, 2, 8, 25)),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 60,vjust = 0.4,hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.y = element_text(hjust=0.5),
        axis.text=element_text(),
        axis.line = element_line(color = "black",size = 0.5, linetype = "solid"))

pdf('/home/dvespasiani/tfbs_plots/tss_eqtls_distance_boxplot.pdf',width = 7,height=7)
tss_distance_plot
dev.off()





