#### script to annotate SNPs with chromHMM data  
library(data.table);library(magrittr);library(dplyr);
library(GenomicRanges);library(R.utils)

setDTthreads(10)
### load cells 
setwd('/data/projects/punim0586/dvespasiani/Files/')

output_dir='./PNG/Chromatin_states/SNPs_chromHMM_annotated/new_set/'
simplified_dir='./PNG/Chromatin_states/simplified_set/new_set/'


read_cells=function(x){
  cells=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y) 
      fread(y,sep = ' ',header=T) %>% 
             makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
             as.data.table()
    )
  
  cell_names=as.character(list.files(x,recursive=F,full.names=F))
  cell_names=gsub("\\..*","",cell_names)
  names(cells)=cell_names
  for (i in seq_along(cells)){
    assign(cell_names[i],cells[[i]],.GlobalEnv)}
  
  return(cells)
  
}

cells=read_cells('./Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states')


read_snps=function(x){
  pop=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)fread(y,sep=' ',header=T,stringsAsFactors = T,
                            select = c('seqnames','start','end','ref','alt',
                                       'ANC','all_frequency','freq_range')) %>% 
             makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
             as.data.table() %>% unique() %>% 
             setnames(c('seqnames','start','end','width','strand','ref','alt','ancestral','all_freq','freq_range'))
    ) 
  pop_names=as.character(list.files(x,recursive=F,full.names=F))
  pop_names=gsub("\\..*","",pop_names)
  names(pop)=pop_names
  for (i in seq_along(pop)){
    assign(pop_names[i],pop[[i]],.GlobalEnv)}
  return(pop)
}

snps=read_snps('./PNG/OoA_snps')

### merge snps with chromatin states ###
lapply(snps,function(x)setkey(x, seqnames, start, end))
lapply(cells,function(x)setkey(x, seqnames, start, end))


chromatin_state_annotation=function(cell,alleles){
  cell=lapply(cell,function(cell)
    alleles=lapply(alleles,function(alleles)
      foverlaps(cell,alleles, type="within")[
        ,element_seqnames :=seqnames
      ][
        ,c('seqnames','i.start','i.end','ref','alt','ancestral','all_freq','freq_range',
           'element_seqnames','start','end','width','chrom_state','cell_type')
           ] %>% 
        setnames(c('seqnames','start','end','ref','alt','ancestral','all_freq','freq_range',
                   'element_seqnames','element_start','element_end','element_width','chrom_state','cell_type'))
    )
  )
}

snps_chromatin_states=chromatin_state_annotation(snps,cells) 

denisova_chromatin_states=rbindlist(snps_chromatin_states[[1]])
neandertal_chromatin_states=rbindlist(snps_chromatin_states[[2]])
nonarchaics_chromatin_states=rbindlist(snps_chromatin_states[[3]])

cell_lines=function(x){x=x[,cell_line := plyr::revalue(cell_type,c("E017"="IMR90_fetal_lung_fibroblast", 
                                                                             
                                                                             "E002"="ES_cells","E008"="ES_cells","E001"="ES_cells",'E015'="ES_cells",'E014'="ES_cells",
                                                                             "E016"="ES_cells", "E003"='ES_cells',"E024"="ES_cells",
                                                                             
                                                                             "E020"="iPSC","E019"="iPSC","E018"="iPSC","E021"="iPSC","E022"="iPSC",
                                                                             
                                                                             "E007"="ES_derived_cells",
                                                                             "E009"="ES_derived_cells","E010"="ES_derived_cells","E013"="ES_derived_cells",
                                                                             "E012"="ES_derived_cells","E011"="ES_derived_cells","E004"="ES_derived_cells",
                                                                             "E005"="ES_derived_cells","E006"="ES_derived_cells",
                                                                             
                                                                             "E062"="TCells","E034"="TCells",
                                                                             "E045"="TCells","E033"="TCells","E044"="TCells","E043"="TCells",
                                                                             "E039"="TCells","E041"="TCells","E042"="TCells","E040"="TCells",
                                                                             "E037"="TCells","E048"="TCells","E038"="TCells","E047"="TCells",
                                                                             
                                                                             "E029"="BCells","E031"="BCells",
                                                                             "E035"="BCells","E051"="BCells","E050"="BCells","E036"="BCells",
                                                                             "E032"="BCells","E046"="BCells","E030"="BCells",
                                                                             
                                                                             "E026"="Mesenchymal","E049"="Mesenchymal",
                                                                             "E025"="Mesenchymal","E023"="Mesenchymal",
                                                                             
                                                                             "E052"='Myosatellite',
                                                                             
                                                                             "E055"="Epithelial","E056"="Epithelial","E059"="Epithelial",
                                                                             "E061"="Epithelial","E057"="Epithelial",
                                                                             "E058"="Epithelial","E028"="Epithelial","E027"="Epithelial",
                                                                             
                                                                             "E054"="Neurospheres","E053"="Neurospheres",
                                                                             
                                                                             "E112"="Thymus",'E093'="Thymus",
                                                                             
                                                                             "E071"="Brain","E074"="Brain",
                                                                             "E068"="Brain","E069"="Brain","E072"="Brain",
                                                                             "E067"="Brain","E073"="Brain","E070"="Brain",
                                                                             "E082"="Brain","E081"="Brain",
                                                                             
                                                                             "E063"="Adipose",
                                                                             
                                                                             "E100"="Muscle","E108"="Muscle","E107"="Muscle","E089"="Muscle","E090"="Muscle",
                                                                             
                                                                             "E083"="Heart","E104"="Heart","E095"="Heart","E105"="Heart","E065"="Heart",
                                                                             
                                                                             "E078"="Smooth_muscle","E076"="Smooth_muscle","E103"="Smooth_muscle","E111"="Smooth_muscle",
                                                                             
                                                                             "E092"="Digestive","E085"="Digestive","E084"="Digestive","E109"="Digestive",
                                                                             "E106"="Digestive","E075"="Digestive","E101"="Digestive","E102"="Digestive",
                                                                             "E110"="Digestive","E077"="Digestive","E079"="Digestive","E094"="Digestive",
                                                                             
                                                                             "E099"="Other_cells","E086"="Other_cells","E088"="Other_cells","E097"="Other_cells",
                                                                             "E087"="Other_cells","E080"="Other_cells",'E091'="Other_cells","E066"="Other_cells",
                                                                             "E098"="Other_cells", "E096"="Other_cells","E113"="Other_cells"
))]
}

snps_chromatin_states=list(denisova_chromatin_states,neandertal_chromatin_states,nonarchaics_chromatin_states) %>% 
  lapply(function(x)cell_lines(x))

snpdensities=function(df,freq){
  if(freq=='nosplit'){
    x=copy(df)
    x=x[
      ,'numbsnps_perelement_perepigenome':= .N, by=.(cell_type,chrom_state,element_seqnames,element_start,element_end)
      ][
        ,'density_snps_perelement_perepigenome':= numbsnps_perelement_perepigenome/element_width]
    
    y=copy(x)
    
    y=y[
      ,c('seqnames','start','end','chrom_state','cell_type','cell_line',
         'density_snps_perelement_perepigenome',
         'numbsnps_perelement_perepigenome')
      ] %>% unique()
    y=y[
      ,'total_snps_perstate_perepigenome':= .N, by=.(cell_type,chrom_state,cell_line)
      ][
        ,'total_snpdens_perepigenome':=sum(density_snps_perelement_perepigenome),by=.(cell_type)
        ][
          ,'sum_snpdens_perstate_perepigenome':=sum(density_snps_perelement_perepigenome),by=.(cell_type,chrom_state)
          ][
            ,'fraction_snpdens_perstate_perepigenome':=sum_snpdens_perstate_perepigenome/total_snpdens_perepigenome]
    z=merge(x,y,all=T)
    return(z)
  } else{
    x=copy(df)
    x=x[
      ,'numbsnps_perelement_perepigenome_perfreqrange':= .N, by=.(cell_type,element_seqnames,element_start,element_end,freq_range)][
        ,'density_snps_perelement_perepigenome_perfreqrange':= numbsnps_perelement_perepigenome_perfreqrange/element_width]
    y=copy(x)
    y=y[
      ,c('seqnames','start','end','chrom_state','cell_type','cell_line',
         'density_snps_perelement_perepigenome_perfreqrange','numbsnps_perelement_perepigenome_perfreqrange','freq_range')
      ] %>% unique()
    y=y[
      ,'total_snps_perstate_perepigenome_perfreqrange':= .N, by=.(cell_type,chrom_state,cell_line,freq_range)
      ][
        ,'total_snpdens_perepigenome_perfreqrange':=sum(density_snps_perelement_perepigenome_perfreqrange),by=.(freq_range,cell_type)
        ][
          ,'sum_snpdens_perstate_perepigenome_perfreqrange':=sum(density_snps_perelement_perepigenome_perfreqrange),by=.(freq_range,cell_type,chrom_state)
          ][
            ,'fraction_snpdens_perstate_perepigenome_perfreqrange':=sum_snpdens_perstate_perepigenome_perfreqrange/total_snpdens_perepigenome_perfreqrange]
    a=merge(x,y,all=T)
    return(a)
  }
}


snps_nofreqsplit=lapply(snps_chromatin_states,function(x)snpdensities(x,'nosplit'))
snps_freqsplit=lapply(snps_chromatin_states,function(x)snpdensities(x,'split'))


combine_data=function(x,y){
  keycolumns=c(
    'seqnames','start','end','ref','alt','ancestral',
    'element_seqnames','element_start','element_end','element_width',
    'chrom_state','cell_type','cell_line','all_freq','freq_range')
  z=purrr::map2(x,y,merge,by=keycolumns,all=T)
  
  pop_names=c('denisova','neandertal','png')
  names(z)=pop_names
  for (i in seq_along(z)){
    assign(pop_names[i],z[[i]],.GlobalEnv)}
    
    return(z)
}

snps_final=combine_data(snps_nofreqsplit,snps_freqsplit) 

filenames=paste0(output_dir,names(snps_final),sep='')
mapply( write.table,snps_final, file = filenames,col.names = T, row.names = F, sep = " ", quote = F)

# simplify dataframes for easy enrichment
simplified=function(x){
  df=copy(x)
  df=df[,.(chrom_state,cell_type,cell_line,freq_range,
       total_snps_perstate_perepigenome,
       fraction_snpdens_perstate_perepigenome,
       total_snps_perstate_perepigenome_perfreqrange,
       fraction_snpdens_perstate_perepigenome_perfreqrange)] %>% unique()
  return(df)
}
snps_final_simplified=lapply(snps_final,function(x)simplified(x))


filenames_simplified=paste0(simplified_dir,names(snps_final_simplified),sep='')
mapply(write.table,snps_final, file = filenames_simplified,col.names = T, row.names = F, sep = " ", quote = F)

