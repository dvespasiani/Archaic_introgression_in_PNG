######################################################
## Table of aSNPs and naSNPs for PNG only deni/nean ##
######################################################
library(dplyr);library(magrittr);library(data.table)
library(ggplot2)
library(openxlsx)

numb_threads=getDTthreads()
threads=setDTthreads(numb_threads-1)

outputdir='./Grouped_filtered_snps/new_set/'
plot_dir='./Results/Plots/QCs/snps_qcs/'
table_output='./Results/Tables/'
setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

##----------------------------------------------------------------------------------
## read files, and for each snp look whether it has a call in archaic haps 
##----------------------------------------------------------------------------------

read_files=function(x){
  x=as.character(list.files(x,recursive =F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep = '\t',header = T)[
        POP_ARCH_REF + POP_ARCH_ALT>0,ancestry := 'archaic' 
        ][
          is.na(ancestry), ancestry := 'non_archaic'
          ])
  
  pop_names=c('deni_png','nean_png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  
}

snps=read_files('./Original_files') 

## look at initial data
nasnps_initial=copy(snps)%>% lapply(function(x)x=x[ancestry=='non_archaic'][!ANC=='-1'])
asnps_initial=copy(snps)%>% lapply(function(x)x=x[ancestry=='archaic'][!ANC=='-1'])

## calculate DAF/MIAF
# calculate_freq=function(x){
#   x=x[
#     ancestry=='archaic',all_frequency := ifelse(POP_ARCH_REF>POP_ARCH_ALT,
#                                                 POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
#                                                 `POP_ARCH_ALT`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
#     ][ancestry=='archaic',state_allele:=ifelse((ANC==0 & POP_ARCH_REF>POP_ARCH_ALT)|(ANC==1 & POP_ARCH_ALT>POP_ARCH_REF),'ancestral',
#                                                ifelse(POP_ARCH_REF==POP_ARCH_ALT,'not_determined','derived'))
#       ][
#         is.na(all_frequency), all_frequency:= ifelse(`ANC`=='0',`POP_NOTARCH_ALT`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`),
#                                                      `POP_NOTARCH_REF`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
#         ][ is.na(state_allele), state_allele:='derived'][
#           ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
#           ]
#   
# }

calculate_freq=function(x){
  x=x[
    ancestry=='archaic',all_frequency := ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                                POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
                                                `POP_ARCH_ALT`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
    ][ancestry=='archaic',state_allele:=ifelse((ANC==0 & POP_ARCH_REF>POP_ARCH_ALT)|(ANC==1 & POP_ARCH_ALT>POP_ARCH_REF),'ancestral',
                                               ifelse(POP_ARCH_REF==POP_ARCH_ALT,'not_determined','derived'))
      ][
        is.na(all_frequency), all_frequency := ifelse(POP_NOTARCH_REF>POP_NOTARCH_ALT,
                                                      POP_NOTARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
                                                      POP_NOTARCH_ALT/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
        ][ is.na(state_allele), state_allele:=ifelse((ANC==0 & POP_NOTARCH_REF>POP_NOTARCH_ALT)|(ANC==1 & POP_NOTARCH_ALT>POP_NOTARCH_REF),'ancestral',
                                                     ifelse(POP_NOTARCH_REF==POP_NOTARCH_ALT,'not_determined','derived'))
           ][
             ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
             ]

}

snps_initial_freq=Map(rbind,asnps_initial,nasnps_initial) %>% lapply(function(x)x=calculate_freq(x))


nasnps_initial_combined=copy(snps_initial_freq)%>% lapply(function(x)x=x[ancestry=='non_archaic']) %>% rbindlist() %>% unique()
  # group_by(`CHR`,`FROM`,`TO`) %>% 
  # filter(`all_frequency`==max(all_frequency)) %>%as.data.table() %>%unique()

all_asnps_initial=copy(snps_initial_freq)%>% lapply(function(x)x=x[ancestry=='archaic']) %>%rbindlist()

nasnps_initial_out_arch_haps=nasnps_initial_combined[!all_asnps_initial,on=c('CHR','FROM','TO')]
nasnps_initial_out_arch_haps=nasnps_initial_out_arch_haps[,pop:='PNG']

asnps_initial=copy(snps_initial_freq) %>% lapply(function(x)x=x[ancestry=='archaic'])
deni_initial=copy(asnps_initial[[1]])[,pop:='Denisova']
nean_initial=copy(asnps_initial[[2]])[,pop:='Neanderthal']

png_initial_list=list(deni_initial,nean_initial,nasnps_initial_out_arch_haps)

##-----------------------
##    plot initial sfs
##-----------------------

get_sfs=function(x){
  df=copy(x)
  df=lapply(df,function(x)x=x[,all_frequency:=round(all_frequency,1)][
    ,numbsnps_perstate_allfreq:=.N,by=.(all_frequency,state_allele)
    ][
      ,totsnps_per_allfreq:=.N,by=.(all_frequency)
      ][,fraction:=round(numbsnps_perstate_allfreq/totsnps_per_allfreq,2)
        ][,log10_numbsnps:=log10(totsnps_per_allfreq)
          ][,c('log10_numbsnps','all_frequency','state_allele','pop','fraction')] %>% unique()) 
 
  return(df) 
}

sfs_initial=get_sfs(png_initial_list)

create_plots=function(x,pop_palette){
  x1=copy(x)[,c('state_allele','fraction'):=NULL] %>% unique()
  x2=copy(x)[,'log10_numbsnps':=NULL] %>% unique()
  
  sfs_palette=c('darkorange1', # ancestral
                        'darkorchid1', #derived
                        'darkolivegreen3' #nd
  )
  
  names(sfs_palette)= levels(as.factor(x$state_allele))
  fraction_col=scale_fill_manual(name= " ",values = sfs_palette,labels = c('Ancestral','Derived','Not determined'))
  
  
  my_palette_pop=c('#C99E10', # denisova
                   '#9B4F0F', #neandertal
                   '#1E656D' #png
  )
  
  names(my_palette_pop)= levels(as.factor(x$pop))
  barcol=scale_fill_manual(name= " ",values = pop_palette)
  
  
  p1=ggplot(x1,aes(x=all_frequency,y=log10_numbsnps,fill=pop))+
    geom_bar(stat = 'identity',position =position_dodge2(width = 0.1))+
    barcol
  
  p2=ggplot(x2,aes(x=all_frequency,y=fraction,fill=state_allele))+
    geom_bar(stat = 'identity',position ='stack')+
    fraction_col
  
  p3=ggpubr::ggarrange(p1,p2,ncol=1)  
return(p3)
}

pdf(paste(plot_dir,'Denisova_initial_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_initial[[1]],'#C99E10')
dev.off()

pdf(paste(plot_dir,'Neanderthal_initial_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_initial[[2]],'#9B4F0F')
dev.off()

pdf(paste(plot_dir,'PNG_initial_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_initial[[3]],'#1E656D')
dev.off()

## plot initial distribution of aSNPs across haplotypes 
distribution_across_haps=function(x){
  x=x[ancestry=='archaic'][,snp_distribution_across_haps :=ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                        (POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT))- (POP_NOTARCH_REF/(POP_NOTARCH_REF+POP_NOTARCH_ALT)),
                                        (POP_ARCH_ALT/(POP_ARCH_REF+POP_ARCH_ALT))-(POP_NOTARCH_ALT/(POP_NOTARCH_REF+POP_NOTARCH_ALT)))
  ][
    ,overall_allele_frequency := ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                        ((POP_ARCH_REF+POP_NOTARCH_REF)/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)),
                                        ((POP_ARCH_ALT+POP_NOTARCH_ALT)/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT)))
    ]
}


deni_initial_distribution=copy(snps_initial_freq[[1]]) %>% distribution_across_haps()
deni_initial_distribution=deni_initial_distribution[,pop:='Denisova']
nean_initial_distribution=copy(snps_initial_freq[[2]]) %>% distribution_across_haps()
nean_initial_distribution=nean_initial_distribution[,pop:='Neanderthal']

asnps_initial_distribution=rbind(deni_initial_distribution,nean_initial_distribution)

pdf(paste(plot_dir,'initial_distribution_across_haps.pdf',sep=''),width = 7,height = 7)
ggplot(asnps_initial_distribution,aes(x=snp_distribution_across_haps,col=pop))+
  geom_density()+
  geom_vline(xintercept=0.25,linetype="dashed")
dev.off()


##-------------------------------------
###       read 1kg frequencies
##-------------------------------------
# read_1kg_snps=function(x){
#   snps=bedr::read.vcf(x)
#   snps_vcf=snps$vcf
#   snps_vcf=snps_vcf[!CHROM%in%c('X','Y')
#            ][,FROM:=POS
#              ][,CHR:= sub("^", "chr", CHROM)
#                ][,paste('1k_pop',1:19,sep='_'):=tstrsplit(INFO,';')
#                  ]
# }
# 
# snps_1kg=read_1kg_snps('../Annotation_and_other_files/human_genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz')
# 
# 
# snps_1kg_modified=copy(snps_1kg)
# snps_1kg_modified=snps_1kg_modified[!CHROM%in%c('X','Y')
#                                     ][,CHR:= sub("^", "chr", CHROM)
#                                       ][,c('CHR','POS','REF','ALT','1k_pop_8','1k_pop_12')] %>%
#   setnames(old=c('POS','1k_pop_8','1k_pop_12'),new = c('FROM','AFR_AF','VT'))
# 
# snps_1kg_modified=snps_1kg_modified[,VT:=gsub(".*=","",VT)][VT%in%'SNP'][,VT:=NULL][,AFR_AF:=gsub(".*=","",AFR_AF)][
#   ,paste('ALT',1:3,sep='_'):=tstrsplit(ALT,',')
#   ][
#     ,paste('AFR_AF',1:3,sep='_'):=tstrsplit(AFR_AF,',')
#     ]
# 
# 
# snps_1kg_multi_alt=copy(snps_1kg_modified) %>% na.omit()
# snps_1kg_alt1=copy(snps_1kg_multi_alt)[,c("CHR",'FROM','REF','ALT_1','AFR_AF_1')] %>% setnames(old=c('ALT_1','AFR_AF_1'),new=c('ALT','AFR_AF'))
# snps_1kg_alt2=copy(snps_1kg_multi_alt)[,c("CHR",'FROM','REF','ALT_2','AFR_AF_2')] %>% setnames(old=c('ALT_2','AFR_AF_2'),new=c('ALT','AFR_AF'))
# snps_1kg_alt3=copy(snps_1kg_multi_alt)[,c("CHR",'FROM','REF','ALT_3','AFR_AF_3')] %>% setnames(old=c('ALT_3','AFR_AF_3'),new=c('ALT','AFR_AF'))
# 
# multi_allele=rbind(snps_1kg_alt1,snps_1kg_alt2,snps_1kg_alt3)
# 
# snps_1kg_no_mult_alt=copy(snps_1kg_modified)
# snps_1kg_no_mult_alt=snps_1kg_no_mult_alt[!multi_allele,on=c('CHR','FROM','REF')]
# snps_1kg_no_mult_alt=snps_1kg_no_mult_alt[,c('CHR','FROM','REF','ALT','AFR_AF')]
# 
# snps_1kg_final=rbind(snps_1kg_no_mult_alt,multi_allele)

# write.table(snps_1kg_final,'../Annotation_and_other_files/human_genome/1kg_snp_modified',sep=' ',quote = F,row.names = F,col.names = T)

snps_1kg_final=fread('../Annotation_and_other_files/human_genome/1kg_snp_modified.gz',sep=' ',header=T)

##-------------------------------------
##  remove SNPs segregating in africa
##-------------------------------------
ooasnps=lapply(snps,function(x)
  x=x[snps_1kg_final,on=c('CHR','FROM','REF','ALT'),nomatch=0][AFR_AF<=0.005][,AFR_AF:=NULL])

snps_notin_1kg=lapply(snps,function(x)x=x[!snps_1kg_final,on=c('CHR','FROM','REF','ALT')])

snps_filtered=Map(rbind,ooasnps,snps_notin_1kg)

##---------------------------------
##     remove singletons 
##---------------------------------
nasnps_no_singletons=copy(snps_filtered) %>% 
  lapply(function(x)x=x[ancestry=='non_archaic'][!ANC=='-1'][
  ,singletons:=ifelse((ANC==0 & POP_NOTARCH_ALT==1)|(ANC==1 & POP_NOTARCH_REF==1),'yes','no') ## removes derived na singletons
                      ][
                        singletons=='no' 
                        ][,singletons:=NULL])

asnps_no_singletons=copy(snps_filtered) %>% 
  lapply(function(x)x=x[ancestry=='archaic'][!ANC=='-1'][
    ,singletons:=
      ifelse((POP_ARCH_REF==0 & POP_ARCH_ALT==1) |(POP_ARCH_REF==1 &POP_ARCH_ALT==0)|(POP_ARCH_REF==1 &POP_ARCH_ALT==1), 'yes','no')
    ][
      singletons=='no' 
      ][,singletons:=NULL])

snps_no_singletons=Map(rbind,asnps_no_singletons,nasnps_no_singletons)

##----------------------------------
##    Calculate MIAF/DAF
##----------------------------------
snp_freq=copy(snps_no_singletons)
snp_freq=lapply(snp_freq,function(x)x=calculate_freq(x))

##----------------------------------------------------------------------------------
##         nasnps are only those consistently outside arch haps
##----------------------------------------------------------------------------------
nasnps_combined=copy(snp_freq)%>% 
  lapply(function(x)x=x[ancestry=='non_archaic']) 
nasnps_outside_arch_haps=nasnps_combined[[1]][nasnps_combined[[2]],on=c('CHR','FROM','TO','REF','ALT','ANC','POP_ARCH_REF','POP_ARCH_ALT',
                                                                        'ancestry','state_allele','DENI_REF','DENI_ALT','NEAN_REF','NEAN_ALT','POP_NOTARCH_REF',
                                                                        'POP_NOTARCH_ALT','all_frequency','freq_range'),nomatch=0]
nasnps=copy(nasnps_outside_arch_haps)
nasnps=nasnps[!(ANC==0 & POP_NOTARCH_ALT==0)][!(ANC==1 & POP_NOTARCH_REF==0)] ## remove nasnps if fixed for ancestral allele

asnps=copy(snp_freq) %>% 
  lapply(function(x)x=x[ancestry=='archaic'][!(POP_ARCH_REF==POP_ARCH_ALT)]) ## remove those for which u cant determine the MIA

##-----------------------------------------------------------------------------------------------------------
##         get ambiguous aSNPs and assign deni/nean when possible assign them to one genome
##-----------------------------------------------------------------------------------------------------------
##                      by looking at archaic reference genome state for each allele 
##------------------------------------------------------------------------------------------------------------

deni_asnps=copy(asnps[[1]])
nean_asnps=copy(asnps[[2]])

ambiguous_snp=inner_join(deni_asnps,nean_asnps,by=c('CHR','FROM','TO','ANC','REF','ALT',"DENI_REF",'DENI_ALT','NEAN_REF','NEAN_ALT')) %>% as.data.table()
ambiguous_snp=ambiguous_snp[
    , all_frequency := ifelse(all_frequency.x > all_frequency.y, all_frequency.x,all_frequency.y)][
      ,POP_ARCH_REF:=ifelse(all_frequency.x >all_frequency.y, POP_ARCH_REF.x,POP_ARCH_REF.y)][
        , POP_ARCH_ALT:=ifelse(all_frequency.x >all_frequency.y,POP_ARCH_ALT.x,POP_ARCH_ALT.y)][
          ,  POP_NOTARCH_REF:=ifelse(all_frequency.x >all_frequency.y,POP_NOTARCH_REF.x,POP_NOTARCH_REF.y)][
            , POP_NOTARCH_ALT:=ifelse(all_frequency.x >all_frequency.y, POP_NOTARCH_ALT.x,POP_NOTARCH_ALT.y)][
              , freq_range :=ifelse(`all_frequency`<0.05,'low','high')
              ] %>% dplyr::select(-c(contains('.y'),contains('.x'))) %>% as.data.table()

ambiguous_snp=ambiguous_snp[
    ,major_introgr_snp := ifelse(POP_ARCH_REF>POP_ARCH_ALT,'ref','alt')][
      ,deni_state :=ifelse(DENI_REF>DENI_ALT,'ref','alt')
      ][
        ,nean_state :=ifelse(NEAN_REF>NEAN_ALT,'ref','alt')
        ][
          ,match_genome := ifelse((major_introgr_snp==deni_state & major_introgr_snp==nean_state) | (major_introgr_snp!=deni_state & major_introgr_snp!=nean_state),
                                  'ambiguous',ifelse((major_introgr_snp==deni_state & major_introgr_snp!=nean_state),'denisova','neandertal'))
          ][
            ,c('major_introgr_snp','deni_state','nean_state'):=NULL
            ]%>% unique()


ambiguous_match_denisova=copy(ambiguous_snp)[match_genome%in%'denisova'][,match_genome:=NULL]
ambiguous_match_neandertal=copy(ambiguous_snp)[match_genome%in%'neandertal'][,match_genome:=NULL]

##---------------------------------------------------------
##         Denisovan/Neanderthal specific aSNPs
##---------------------------------------------------------
archaic_specific=function(x,y){
  df=copy(x)
  df=df[!ambiguous_snp, on=c("CHR", "FROM",'TO')]
  y=y[,ancestry:='archaic'][,state_allele:=ifelse((ANC==0 & POP_ARCH_REF>POP_ARCH_ALT)|(ANC==1 & POP_ARCH_ALT>POP_ARCH_REF),'ancestral',
                                            ifelse(POP_ARCH_REF==POP_ARCH_ALT,'not_determined','derived'))]
  df=rbind(df,y) %>% unique()
  return(df)
  
}

denisova_asnps_plus_ambiguous=archaic_specific(asnps[[1]],ambiguous_match_denisova)
neandertal_asnps_plus_ambiguous=archaic_specific(asnps[[2]],ambiguous_match_neandertal)

## remove aSNPs that are broadly distributed across haplotypes
## set threshold for aSNPs to >=0.25, i.e. aSNPs must occur more than 25% of the times within arch haps 
denisova_asnps=distribution_across_haps(denisova_asnps_plus_ambiguous)[,pop:='Denisova']
neandertal_asnps=distribution_across_haps(neandertal_asnps_plus_ambiguous)[,pop:='Neanderthal']

pdf(paste(plot_dir,'filtered_snps_distribution_across_haps.pdf',sep=''),width = 7,height = 7)
ggplot(rbind(denisova_asnps,neandertal_asnps),aes(x=snp_distribution_across_haps,col=pop))+
  geom_density()+
  geom_vline(xintercept=0.25,linetype="dashed")
dev.off()


denisova_asnps_arch_haps=copy(denisova_asnps)[
  ,snp_distribution_across_haps:=round(snp_distribution_across_haps,2)
  ][snp_distribution_across_haps>=0.25][,c('snp_distribution_across_haps','overall_allele_frequency'):=NULL]

neandertal_asnps_arch_haps=copy(neandertal_asnps)[
  ,snp_distribution_across_haps:=round(snp_distribution_across_haps,2)
  ][snp_distribution_across_haps>=0.25][,c('snp_distribution_across_haps','overall_allele_frequency'):=NULL]


png_final_list=list(denisova_asnps_arch_haps,neandertal_asnps_arch_haps,nasnps)
names(png_final_list)=c('denisova','neanderthal','png')

##-----------------------------------------------------------------------------------------------------------
##  before writing the files remember to get the frequency of the alternative allele in the png naSNP set
##-----------------------------------------------------------------------------------------------------------

nasnps_daf=copy(nasnps)
nasnps_daf=nasnps_daf[
  ,all_frequency:= ifelse(`ANC`=='0',
                          POP_NOTARCH_ALT/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`),
                           `POP_NOTARCH_REF`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
  ][
    ,state_allele:='derived'
    ][
      ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
      ]

png_list_towrite=list(denisova_asnps_arch_haps,neandertal_asnps_arch_haps,nasnps_daf)
names(png_list_towrite)=c('denisova','neanderthal','png')

# filenames_png_snps=paste(outputdir,names(png_list_towrite),sep='')
# mapply(write.table,png_list_towrite, file = filenames_png_snps,col.names = T, row.names = F, sep = " ", quote = F)

##-----------------------
##    plot final sfs
##-----------------------
png_final_list=Map(mutate,png_final_list,'pop'=names(png_final_list))
png_final_list=png_final_list %>% lapply(function(x)x=x %>% as.data.table(x))


sfs_final=get_sfs(png_final_list)

pdf(paste(plot_dir,'Denisova_final_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_final[[1]],'#C99E10')
dev.off()

pdf(paste(plot_dir,'Neanderthal_final_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_final[[2]],'#9B4F0F')
dev.off()

pdf(paste(plot_dir,'PNG_final_sfs.pdf',sep=''),width = 15,height = 7)
create_plots(sfs_final[[3]],'#1E656D')
dev.off()


##---------------------------------------------------------
##        Create Supplementary Table
##---------------------------------------------------------
count_snps=function(x){x=x[,c('CHR','FROM','TO')] %>% unique() %>% nrow()}

make_table=function(x){
  
  snp_per_file=copy(x)
  
  common=copy(x)
  common=common[[1]][common[[2]],on=c('CHR','FROM','TO','REF','ALT'),nomatch=0] %>% count_snps()
  
  snp_per_file=lapply(snp_per_file,function(y)count_snps(y))
  
  dt=data.table('Denisova file'=snp_per_file[[1]], 'Neanderthal file'=snp_per_file[[2]],'Common SNPs'=common)
    
  return(dt)
}

common_column=function(x){
  common=copy(x)
  common=common[[1]][common[[2]],on=c('CHR','FROM','TO','REF','ALT'),nomatch=0]
  
}


final_table=function(...){
  x=list(...) %>% rbindlist()
  
  objects=as.list(substitute(list(...)))[-1L]
  objects=setNames(list(...), objects)
  
  x=x[,'Variable':=names(objects)] %>% dplyr::select(c('Variable',everything())) %>% as.data.table()
  return(x)
}

## first sheet
total_snps=make_table(snps)

tot_asnps=copy(snps) %>% lapply(function(x)x=x[ancestry=='archaic'])
tot_asnps=make_table(tot_asnps)

tot_nasnps=copy(snps) %>% lapply(function(x)x=x[ancestry=='non_archaic'])
tot_nasnps=make_table(tot_nasnps)

tot_asnps_ooa=copy(snps_filtered)%>% lapply(function(x)x=x[ancestry=='archaic'])
tot_asnps_ooa=make_table(tot_asnps_ooa)

tot_nasnps_ooa=copy(snps_filtered)%>% lapply(function(x)x=x[ancestry=='non_archaic'])
tot_nasnps_ooa=make_table(tot_nasnps_ooa)

tot_asnps_no_anc_info=copy(snps_filtered) %>% lapply(function(x)x=x[ancestry=='archaic'& ANC==-1])
tot_asnps_no_anc_info=make_table(tot_asnps_no_anc_info)

tot_nasnps_no_anc_info=copy(snps_filtered) %>% lapply(function(x)x=x[ancestry=='non_archaic'& ANC==-1])
tot_nasnps_no_anc_info=make_table(tot_nasnps_no_anc_info)

tot_asnps_no_singletons=copy(snps_filtered)%>% lapply(function(x)x=x[ancestry=='archaic'][!((POP_ARCH_REF==0 & POP_ARCH_ALT==1) |(POP_ARCH_REF==1 &POP_ARCH_ALT==0)|(POP_ARCH_REF==1 &POP_ARCH_ALT==1))])
tot_asnps_no_singletons=make_table(tot_asnps_no_singletons)

tot_nasnps_no_singletons=copy(snps_filtered)%>% lapply(function(x)x=x[ancestry=='non_archaic'][(ANC==0 & POP_NOTARCH_ALT==1)|(ANC==1 & POP_NOTARCH_REF==1)])
tot_nasnps_no_singletons=make_table(tot_nasnps_no_singletons)

tot_asnps_miaf_assigned=copy(snp_freq) %>%lapply(function(x)x=x[ancestry=='archaic'][!(POP_ARCH_REF==POP_ARCH_ALT)])
tot_asnps_miaf_assigned=make_table(tot_asnps_miaf_assigned)

tot_nasnps_non_fixed_ancestral=copy(snp_freq) %>%lapply(function(x)x=x[ancestry=='non_archaic'][!((ANC==0 & POP_NOTARCH_ALT==0)|(ANC==1 & POP_NOTARCH_REF==0))])
tot_nasnps_non_fixed_ancestral=make_table(tot_nasnps_non_fixed_ancestral)


first_sheet=final_table(total_snps,tot_asnps,tot_nasnps,tot_asnps_ooa,tot_nasnps_ooa,
                        tot_asnps_no_anc_info,tot_nasnps_no_anc_info,
                        tot_asnps_no_singletons,tot_nasnps_no_singletons,
                        tot_asnps_miaf_assigned,tot_nasnps_non_fixed_ancestral)
                          
## second sheet
second_make_table=function(x){
 number_snps=copy(x)
 number_snps=count_snps(number_snps)
 dt=data.table('Number of SNPs'=number_snps)
 return(dt)
}


second_final_table=function(...){
  x=list(...) %>% rbindlist()
  
  objects=as.list(substitute(list(...)))[-1L]
  objects=setNames(list(...), objects)
  
  x=x[,'Variable':=names(objects)] %>% dplyr::select(c('Variable',everything())) %>% as.data.table()
  return(x)
}

tot_deni_plus_ambig=second_make_table(denisova_asnps)
tot_nean_plus_ambig=second_make_table(neandertal_asnps)

tot_deni_mainly_arch_haps=second_make_table(denisova_asnps_arch_haps)
tot_nean_mainly_arch_haps=second_make_table(neandertal_asnps_arch_haps)
tot_nasnps=second_make_table(nasnps_daf)
  
tot_rare_deni_asnps=copy(denisova_asnps_arch_haps)[freq_range=='low'] %>% second_make_table()
tot_rare_nean_asnps=copy(neandertal_asnps_arch_haps)[freq_range=='low'] %>% second_make_table()
tot_rare_nasnps=copy(nasnps_daf)[freq_range=='low'] %>% second_make_table()


second_sheet=second_final_table(tot_deni_plus_ambig,tot_nean_plus_ambig,
                                tot_deni_mainly_arch_haps,tot_nean_mainly_arch_haps,tot_nasnps,
                                tot_rare_deni_asnps,tot_rare_nean_asnps,tot_rare_nasnps)


supp_table=list(first_sheet,second_sheet)
write.xlsx(supp_table,paste(table_output,'Supp_Table_1.xlsx',sep=''))
