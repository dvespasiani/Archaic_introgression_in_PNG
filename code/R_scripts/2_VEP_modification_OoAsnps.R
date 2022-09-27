### Combine VEP output and remove SNPs shared with Africans 
library(tidyr);library(data.table);library(magrittr);library(dplyr);
library(R.utils)

setDTthreads(8)
setwd('/data/projects/punim0586/dvespasiani/Files/Archaic_introgression_in_PNG/')

output_dir='./OoA_snps/'

read_files=function(x){
  snps=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)y=fread(y,sep='\t',header = T))
  snps=lapply(snps,function(y)y=y[
    ,c('seqnames','start','allele'):= tstrsplit(`#Uploaded_variation` , "_", fixed=TRUE)
    ][
      ,c('ref','alt'):= tstrsplit(`allele` , "/", fixed=TRUE)
      ][
        ,genomic_element:=gsub("\\.,*","",Consequence)
        ][
          ,c('allele','Allele','Consequence','Location','#Uploaded_variation'):=NULL
          ])
  snps=lapply(snps,function(x)x=x[,start:=as.numeric(start)][,end:=start+1])
  return(snps)
}

assign_names=function(x){
  pop_names=c('denisova','neandertal','png')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  return(x)
}

vep_output=read_files('./VEP/output/')
vep_output=assign_names(vep_output)
## add snp frequencies from original files 
read_original_snps=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,header = T,sep = ' ') [
        ,freq_range:= ifelse(`all_frequency`<0.05,'low','high')
        ]
    )
}
original_snps=read_original_snps('./Grouped_filtered_snps/new_set/')

combined_files=purrr::map2(vep_output,original_snps,inner_join,by=c('seqnames'='CHR','start'='FROM','end'='TO','ref'='REF','alt'='ALT')) 
combined_files=Map(mutate,combined_files,'pop'=names(combined_files)) %>% lapply(function(x)setDT(x))

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

ooasnps=out_of_africa_snps(combined_files)
# lapply(ooasnps,function(y)y=y[,c('seqnames','start','end','AF')]%>%unique()%>%nrow())

filenames=paste0(output_dir,names(ooasnps),sep='')
mapply(write.table,ooasnps, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)



