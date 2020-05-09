## D1 enrichment over D2 in PNG ##
library(data.table);library(magrittr);library(dplyr);
library(GenomicRanges);
library(R.utils)

setDTthreads(8)
# READ COMPONENTS 
setwd('/data/projects/punim0586/dvespasiani/Files/') 

read_components=function(x){
  x=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y)
      fread(y,sep='\t',select = c(1:3),col.names = c('chr','start','end'),header = F)[
        ,seqnames:= sub("^", "chr", chr)
        ][
          ,chr:=NULL
          ] )
  pop_names=c('1','2')
  names(x)=pop_names
  for (i in seq_along(x)){
    assign(pop_names[i],x[[i]],.GlobalEnv)}
  x=Map(mutate,x,'deni_component'=names(x)) %>% rbindlist()
  return(x)
}

denisova_components=read_components('./Denisova_components/')

# Denisova SNPs 
denisova_in_png=fread('./PNG/OoA_snps/denisova_png',sep=' ',header = T)[
  ,c('seqnames','start','end','ref','alt','ANC','all_frequency')
  ][
    ,freq_range:=ifelse(all_frequency<0.05,'low','high')
    ][
      ,all_freq:=all_frequency
      ] [
        ,ancestral:=ANC
        ][
          ,c('ANC','all_frequency'):=NULL
          ]%>% unique()

# merge deni snps with components
keys=c('seqnames','start','end')
setkeyv(denisova_in_png,keys)
setkeyv(denisova_components,keys)

snp_component=foverlaps(denisova_in_png,denisova_components,type='within')[
  ,c('start','end'):=NULL
  ][
    ,start:=i.start
    ][
      ,end:=i.end
      ][
        ,c('i.start','i.end'):=NULL
        ]%>% mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
  mutate_at('deni_component',as.factor) %>% 
  as.data.table() %>% unique()


denisovans=snp_component %>% split(as.factor(snp_component$deni_component))
denisovans=lapply(denisovans,function(x)x=x[,'deni_component':=NULL]%>% unique())

# some snps are both in D1 and D2 haps -> assign them to D0
unassigned_denisova=semi_join(denisovans[[2]],denisovans[[3]],by=keys) # 531 are shared between D1 and D2 haps

indo_denisova=anti_join(denisovans[[3]],unassigned_denisova,by=keys)
png_denisova=anti_join(denisovans[[2]],unassigned_denisova,by=keys)

isea_denisovans=rbind(indo_denisova,png_denisova)

altai_denisova=anti_join(denisovans[[1]],isea_denisovans,by=keys)
altai_denisova=rbind(altai_denisova,unassigned_denisova)

sorted_denisovans=list(altai_denisova,png_denisova,indo_denisova) %>% 
  lapply(function(x)as.data.table(x))
# lapply(sorted_denisovans,function(x)
# x=x[,c('seqnames','start','end')] %>% unique() %>% nrow())
# 5274 D1 snps (i.e. ~ 4%)
# 5193 D2 snps (i.e. ~ 4% )
# 133114 D0 snps (i.e. ~ 92%)
# 143581 total Denisova aSNPs

## annotate them again with chromstates
assign_cell_names=function(x,list_df){
  x=as.character(list.files(x,recursive=F,full.names=F))
  names(list_df)=x
  for (i in seq_along(list_df)){
    assign(x[i],list_df[[i]],.GlobalEnv)}
  return(list_df)
}

read_cells=function(x){
  df=as.character(list.files(x,recursive = F,full.names = T)) %>% 
    lapply(function(y) fread(y,sep = ' ',header=T) %>% makeGRangesFromDataFrame(keep.extra.columns =T) %>% 
             as.data.table()
    )
  df=assign_cell_names(x,df)
  return(df)
}


cells=read_cells('./Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states')

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

lapply(sorted_denisovans,function(x)setkey(x, seqnames, start, end))
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
                   'element_seqnames','element_start','element_end','element_width',
                   'chrom_state','cell_type')) %>% unique()
    )
  )
}

snps_chromatin_states=chromatin_state_annotation(sorted_denisovans,cells) 

altai_chromatin_states=rbindlist(snps_chromatin_states[[1]])
d1_chromatin_states=rbindlist(snps_chromatin_states[[2]])
d2_chromatin_states=rbindlist(snps_chromatin_states[[3]])

denisovans_chromatin_states=list(altai_chromatin_states,d1_chromatin_states,d2_chromatin_states)
denisovans_chromatin_states=lapply(denisovans_chromatin_states,function(x)cell_lines(x))

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


snps_nofreqsplit=lapply(denisovans_chromatin_states,function(x)snpdensities(x,'nosplit'))
snps_freqsplit=lapply(denisovans_chromatin_states,function(x)snpdensities(x,'split'))

#merge dfs
combine_data=function(x,y){
  keycolumns=c(
    'seqnames','start','end','ref','alt','ancestral',
    'element_seqnames','element_start','element_end','element_width',
    'chrom_state','cell_type','cell_line','all_freq','freq_range')
  z=purrr::map2(x,y,merge,by=keycolumns,all=T)
  pop_names=c('0','1','2')
  names(z)=pop_names
  for (i in seq_along(z)){
    assign(pop_names[i],z[[i]],.GlobalEnv)}
  
  z=Map(mutate,z,'deni_component'=names(z)) %>% rbindlist()
  return(z)
}

snps_final=combine_data(snps_nofreqsplit,snps_freqsplit) 

deni_output_dir='./PNG/Chromatin_states/D1_D2/new_set/'
filenames=paste0(deni_output_dir,'denisovans_snpdensities',sep='')
simplified_filnames=paste0(deni_output_dir,'denisovans_simplified',sep='')

write.table(snps_final,filenames,row.names=FALSE, sep=" ", quote=F)

## need only chromstate + cell + freq + densities + deni component
simplified=function(x){
  x[,.(chrom_state,cell_type,cell_line,freq_range,deni_component,
       total_snps_perstate_perepigenome,
       fraction_snpdens_perstate_perepigenome,
       total_snps_perstate_perepigenome_perfreqrange,
       fraction_snpdens_perstate_perepigenome_perfreqrange)] %>% unique()
}
snps_simplified=simplified(snps_final)

write.table(snps_simplified,simplified_filnames,row.names=FALSE, sep=" ", quote=F)







