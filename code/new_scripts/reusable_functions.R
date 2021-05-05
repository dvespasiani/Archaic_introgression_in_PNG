## tissue 
tissue = c(
    'Adipose','BCells','Brain',
    'Digestive','Epithelial','ES_cells',
    'ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
    'iPSC','Mesenchymal','Muscle','Myosatellite',
    'Neurospheres','Other_cells','Smooth_muscle',
    'TCells','Thymus'
)

## make tissue column
group_celltypes = function(x){x=x[,cell_line := plyr::revalue(cell_type,c("E017"="IMR90_fetal_lung_fibroblast", 
                                                                   
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

## calculate odds ratio
calculate_odds_ratio = function(asnps,nasnps,merging_keys,fulljoin){
    table = merge(asnps,nasnps,by=c(merging_keys),all=fulljoin)
    table[is.na(table)] = 0
    table = split(table,by=merging_keys)%>%
    lapply(function(x)x=x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%
    fisher.test()
    )
    
    fisher_test_results = copy(table)%>%
    lapply(function(x)x=data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
    ))
    fisher_test_results = Map(mutate,fisher_test_results,elements=names(fisher_test_results))%>%rbindlist()
    fisher_test_results = fisher_test_results[,significance:=ifelse(p<0.05,'*','')]

    return(fisher_test_results)
}

## read tfbs SNPs
read_tfbs = function(dir,patterns){
  tfbs = list.files(dir,recursive = T,full.names = T,pattern = patterns) %>% 
  lapply(
    function(y)fread(y,sep='\t',header = T)[,end:=start+1]
    )
    names(tfbs) = gsub("\\.txt$","",basename(list.files(dir,recursive = T,full.names = F,pattern=patterns)))
  return(tfbs)
}

## read cres snps
read_cres = function(dir,cols,cells){
    cres = list.files(dir,recursive = F,full.names = T) %>%
    lapply(
        function(y) z=fread(y,sep='\t',header = T,select = cols)[
            chrom_state %in% c('1_TssA','2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
            ][
            ,MAF:=round(MAF,2)
            ][
                ,chrom_state:=NULL
                ]%>%unique()%>%setnames(old=7,new = 'allele_of_interest')
        )
    if(missing(cells)) {
        cres_allcells = lapply(cres,function(x)x=x[cell_line %in% tissue])
    return(cres_allcells)
    } else {
        cres_selected_cells = lapply(cres,function(x)x=x[cell_line %in% cells])
    return(cres_selected_cells)
    }
}