### Combine VEP output and remove SNPs shared with Africans 
library(tidyr);library(data.table);library(magrittr);library(dplyr);
library(R.utils)

setDTthreads(8)
setwd('/data/projects/punim0586/dvespasiani/Files/PNG/')

output_dir='./OoA_snps/'

read_asnps_files=function(x){
  asnps=as.character(list.files(x,recursive = F,full.names = T))
  asnps=asnps[c(2,3)] %>%# dont consider ambiguous and separately treat naSNPs because you have already modified the df in another script
    lapply(function(y)y=fread(y,sep='\t',header = T))
  asnps=lapply(asnps,function(y)y=y[
    ,c('seqnames','start','allele'):= tstrsplit(`#Uploaded_variation` , "_", fixed=TRUE)
    ][
      ,c('ref','alt'):= tstrsplit(`allele` , "/", fixed=TRUE)
      ][
        ,c('genomic_element','Consequence') :=tstrsplit(`Consequence` , ",", fixed=TRUE)
        ][
          ,c('allele','Allele','Consequence','Location','#Uploaded_variation'):=NULL
          ][
            ,start := as.numeric(start)
            ][
              ,end:= as.numeric(start)+1
              ] %>% unique())
  
  return(asnps)
}

assign_names=function(x){
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

vep_output=read_asnps_files('./VEP/output/')
nasnps_vep=fread('./VEP/output/non_archaics_png.txt.gz',sep=' ',header=T)

vep_output=list(vep_output[[1]],vep_output[[2]],nasnps_vep) %>% assign_names()

## add snp frequencies from original files 
read_original_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,header = T,sep = ' ') [
        ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
        ]
    )
  x=x[c(2:4)]
}
original_snps=read_original_snps('./Grouped_filtered_snps')

combined_files=purrr::map2(vep_output,original_snps,inner_join,by=c('seqnames'='CHR','start'='FROM','end'='TO','ref'='REF','alt'='ALT')) 
combined_files=Map(mutate,combined_files,'pop'=names(combined_files)) %>% lapply(function(x)setDT(x))
#numb snp ambig=91785;deni=475422;nean=310680;png=4407553


out_of_africa_snps=function(x){
  df=copy(x)
  df=lapply(df,function(y)
    y=y[
      ,AFR_AF :=ifelse(AFR_AF=='-',as.numeric(0),as.numeric(AFR_AF))
      ][
        ,c('overall_allele_frequency','snp_distribution_across_haps'):=NULL
        ][AFR_AF<0.005])
  return(df)
}

ooasnps=out_of_africa_snps(combined_files)## % snps removed is between 20% ambig; 30% deni; 37% nean; 25% png 
# lapply(ooasnps,function(y)y=y[,c('seqnames','start','end')]%>%unique()%>%nrow())

filenames=paste0(output_dir,names(ooasnps),sep='')
mapply(write.table,ooasnps, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)



