######################################################
## Table of aSNPs and naSNPs for PNG only deni/nean ##
######################################################
library(dplyr);library(magrittr);library(data.table)
library(GenomicRanges)
library(ggplot2)

setDTthreads(10)

outputdir='./Grouped_filtered_snps/'

vep_output_dir= './VEP/input/'
  
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

read_files=function(x){
  x=as.character(list.files(x,recursive =F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep = '\t',header = T)[
        POP_ARCH_REF + POP_ARCH_ALT>0,ancestry := 'archaic' 
        ][
          is.na(ancestry), ancestry := 'non_archaic'
          ][
            ancestry=='archaic',all_frequency := ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                                        POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
                                                        `POP_ARCH_ALT`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
            ][
              is.na(all_frequency), all_frequency :=ifelse(`ANC`=='0',`POP_NOTARCH_ALT`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`),
                                                           `POP_NOTARCH_REF`/(`POP_ARCH_REF`+`POP_ARCH_ALT`+`POP_NOTARCH_REF`+`POP_NOTARCH_ALT`))
              ][
                ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
                ][
                  !(ancestry=='non_archaic' & ANC=='-1')]     
    )
  
  pop_names=c('deni_png','nean_png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
  
}


population_files=read_files('./Original_files') 

number_snps=function(x){df=copy(x)[,c(1:3)] %>%unique() %>%  nrow()}

tot_snps_between_haplotypes=copy(population_files) %>% rbindlist() %>% number_snps() # 8192985
tot_snps_per_file_per_haplotypes=lapply(population_files,function(x) x=x %>% split(as.factor(x$ancestry)) %>% lapply(function(y) y%>% number_snps())) 

## separate aSNPs from naSNPs and remove singletons 

remove_singeletons=function(file,population){
  x=copy(file)
  if(population%in%'archaic'){
    x=lapply(x,function(y)y=y[ancestry%in%population][POP_ARCH_REF+POP_ARCH_ALT>1]) # removes all archaic singletons
  }else{
    x=lapply(x,function(y)y=y[ancestry%in%population][ ANC==0 & POP_NOTARCH_ALT >1 |ANC==1 & POP_NOTARCH_REF>1 ])# removes derived singletons
  }
  return(x)
  
}

aSNPs_nosingletons=remove_singeletons(population_files,'archaic')
naSNPs_nosingletons=remove_singeletons(population_files,'non_archaic')

tot_nasnps_common_between_files=population_files[[1]][ancestry%in%'non_archaic'][population_files[[2]][ancestry%in%'non_archaic'],on=c('CHR',"FROM",'TO'),nomatch=0] %>% number_snps() # 5874240
tot_nasnps_nosingletons_common_between_files=naSNPs_nosingletons[[1]][naSNPs_nosingletons[[2]],on=c('CHR',"FROM",'TO'),nomatch=0] %>% number_snps() # 4153590

# combine naSNPs and remove those SNPs that have instances in any of the 2 archaic files
list_of_nasnps=function(x){
  df=copy(x) %>% rbindlist()%>% 
    group_by(`CHR`,`FROM`,`TO`) %>% filter(`all_frequency`==max(all_frequency)) %>%as.data.table() %>%unique()
  
  asnps_combined=copy(aSNPs_nosingletons) %>% rbindlist()
  
  df=df[!asnps_combined, on=c("CHR", "FROM",'TO')]
  return(df)
}

naSNPs_list=list_of_nasnps(naSNPs_nosingletons) 

tot_nasnps_combined=number_snps(naSNPs_list) # 4407553


# get ambiguous + deni/nean specific aSNPs

denisova_aSNPs=aSNPs_nosingletons[[1]]
neandertal_aSNPs=aSNPs_nosingletons[[2]]


ambiguous_variants=function(x,y){
  denisova_copy=copy(x)
  neandertal_copy=copy(y)
  
  denisova_copy=denisova_copy[,c(11:15):=NULL]
  ambiguous=inner_join(denisova_copy,neandertal_copy,by=c('CHR','FROM','TO','ANC','REF','ALT')) %>% 
    as.data.table()
  ambiguous=ambiguous[
    , all_frequency := ifelse(all_frequency.x > all_frequency.y, all_frequency.x,all_frequency.y)][
      ,POP_ARCH_REF:=ifelse(all_frequency.x >all_frequency.y, POP_ARCH_REF.x,POP_ARCH_REF.y)][
        , POP_ARCH_ALT:=ifelse(all_frequency.x >all_frequency.y,POP_ARCH_ALT.x,POP_ARCH_ALT.y)][
          ,  POP_NOTARCH_REF:=ifelse(all_frequency.x >all_frequency.y,POP_NOTARCH_REF.x,POP_NOTARCH_REF.y)][
            , POP_NOTARCH_ALT:=ifelse(all_frequency.x >all_frequency.y, POP_NOTARCH_ALT.x,POP_NOTARCH_ALT.y)][
              , freq_range :=ifelse(`all_frequency`<0.05,'low','high')
              ][
                ,c('CHR','FROM','TO','ANC','REF','ALT',"all_frequency","freq_range","POP_ARCH_REF",
                   "POP_ARCH_ALT","POP_NOTARCH_REF","POP_NOTARCH_ALT",
                   "DENI_REF","DENI_ALT","NEAN_REF","NEAN_ALT","ancestry"       
                )
                ]
  # # look at archaic reference genome state for each allele to try and assign ambiguous to one genome
  ambiguous=ambiguous[
    ,major_introgr_snp := ifelse(POP_ARCH_REF>POP_ARCH_ALT,'ref','alt')][
      ,deni_state :=ifelse(DENI_REF>DENI_ALT,'ref','alt')
      ][
        ,nean_state :=ifelse(NEAN_REF>NEAN_ALT,'ref','alt')
        ][
          ,match_genome := ifelse(major_introgr_snp==deni_state & major_introgr_snp==nean_state | major_introgr_snp!=deni_state & major_introgr_snp!=nean_state,'ambiguous',
                                  ifelse(major_introgr_snp==deni_state & major_introgr_snp!=nean_state,'denisova','neandertal'))
          ][
            ,c('major_introgr_snp','deni_state','nean_state'):=NULL
            ]%>% unique()
  return(ambiguous)
  
}

ambiguous_aSNPs=ambiguous_variants(denisova_aSNPs,neandertal_aSNPs)

ambiguous_aSNPs_specific=copy(ambiguous_aSNPs)[match_genome%in%'ambiguous'][,match_genome:=NULL]
ambiguous_match_denisova=copy(ambiguous_aSNPs)[match_genome%in%'denisova'][,match_genome:=NULL]
ambiguous_match_neandertal=copy(ambiguous_aSNPs)[match_genome%in%'neandertal'][,match_genome:=NULL]


tot_ambiguous_aSNPs=number_snps(ambiguous_aSNPs_specific) # 219403
tot_ambiguous_match_deni=number_snps(ambiguous_match_denisova) # 16231
tot_ambiguous_match_nean=number_snps(ambiguous_match_neandertal) # 19641

# denisova/neandertal specific snps
archaic_specific=function(x,y){
  df=copy(x)
  df=df[!ambiguous_aSNPs, on=c("CHR", "FROM",'TO')]
  
  df=rbind(df,y) %>% unique()
  return(df)
  
}

denisova_aSNPs_specific=archaic_specific(denisova_aSNPs,ambiguous_match_denisova)
neandertal_aSNPs_specific=archaic_specific(neandertal_aSNPs,ambiguous_match_neandertal)

## remove aSNPs that are broadly distributed across haplotypes
## set threshold for aSNPs to 0.25, i.e. aSNPs must occur more than 25% of the times within arch haps 
asnps_between_haps=function(x){
  df=copy(x)
  df=df[,snp_distribution_across_haps :=ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                               POP_ARCH_REF/(POP_ARCH_REF+POP_ARCH_ALT)- POP_NOTARCH_REF/(POP_NOTARCH_REF+POP_NOTARCH_ALT),
                                               POP_ARCH_ALT/(POP_ARCH_REF+POP_ARCH_ALT)- POP_NOTARCH_ALT/(POP_NOTARCH_REF+POP_NOTARCH_ALT))
        ][
          ,overall_allele_frequency := ifelse(POP_ARCH_REF>POP_ARCH_ALT,
                                              (POP_ARCH_REF+POP_NOTARCH_REF)/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
                                              (POP_ARCH_ALT+POP_NOTARCH_ALT)/(POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT))
          ][snp_distribution_across_haps>0.25]
  return(df)
}

denisova_aSNPs_within_arch_haplotypes=asnps_between_haps(denisova_aSNPs_specific)
neandertal_aSNPs_within_arch_haplotypes=asnps_between_haps(neandertal_aSNPs_specific)

tot_deni_snps=number_snps(denisova_aSNPs_within_arch_haplotypes)
tot_nean_snps=number_snps(neandertal_aSNPs_within_arch_haplotypes)


png_snps_list=list(denisova_aSNPs_within_arch_haplotypes,neandertal_aSNPs_within_arch_haplotypes,naSNPs_list)

assign_names=function(x){
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

png_snps_list=assign_names(png_snps_list)

# filenames=paste0(outputdir,names(png_snps_list),sep='')
# mapply(write.table,png_snps_list, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)

## output files for VEP annotation
vep_input_format=copy(png_snps_list) %>% lapply(function(x)x[,c('CHR','FROM','REF','ALT')][
  ,TO := FROM
  ][
    ,strand:= '+'
    ][
      ,allele := paste(REF,ALT,sep='/')
      ][
        ,c('CHR','FROM','TO','strand','allele')
        ] %>% unique())


vep_filenames=paste0(vep_output_dir,names(vep_input_format),sep='')
vep_filenames=paste0(vep_filenames,'.txt',sep='')

# mapply(write.table,vep_input_format, file = vep_filenames,col.names = T, row.names = F, sep = "\t", quote = F)


## plot sfs
sfs=copy(png_snps_list)
sfs=lapply(sfs,function(x)x=x[,all_frequency:=round(all_frequency,1)][
  ,numbsnps:=.N,by=all_frequency
][,log10_numbsnps:=log10(numbsnps)][,c('log10_numbsnps','all_frequency')] %>% unique())

sfs=Map(mutate,sfs,'pop'=names(sfs)) %>% rbindlist()

pdf('/home/dvespasiani/sfs.pdf',width = 7,height = 7)
ggplot(sfs,aes(x=all_frequency,y=log10_numbsnps,fill=pop))+
  geom_bar(stat = 'identity',position = 'dodge')
dev.off()


