##-----------
## vectors 
##-----------
range_keys = c('seqnames','start','end')
chrom_states = c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")
cres_states = chrom_states[c(1:3,6,7)] 
ancestry <- c('denisova','modern_humans','neanderthal')
## color palette
my_palette_pop=c(
    '#82c0cc', # denisova
    '#ede7e3', # modern humans
    '#ffa62b'  # neandertal
)
names(my_palette_pop) = ancestry

## tissue 
tissue = c(
    'Adipose','BCells','Brain',
    'Digestive','Epithelial','ES_cells',
    'ES_derived_cells','Heart','IMR90_fetal_lung_fibroblast',
    'iPSC','Mesenchymal','Muscle','Myosatellite',
    'Neurospheres','Other_cells','Smooth_muscle',
    'TCells','Thymus'
)

tissue_colors = c(
  "Adipose"='darkorange3',
  "BCells"='palegreen4',
  "Brain"='goldenrod3',
  "Digestive"="plum3",
  "Epithelial"="orange1",
  "ES_cells"='hotpink4',
  "ES_derived_cells"='dodgerblue2',
  "Heart"="palevioletred2",
  'IMR90_fetal_lung_fibroblast'='red3',
  "iPSC"='mediumpurple4',
  "Mesenchymal"='hotpink3',
  "Muscle"="indianred",
  "Myosatellite"='darkorange2',
  "Neurospheres"='lightgoldenrod2',
  "Other_cells"="grey60",
  "Smooth_muscle"="hotpink1",
  "TCells"='chartreuse4',
  "Thymus"='yellow2'
)

chromHMM_palette = c(
  '#FF0000','#FF6E00','#32CD32',
  '#008000','#006400','#C2E105',
  '#FFFF00','#66CDAA','#8A91D0',
  '#CD5C5C','#E9967A','#BDB76B',
  '#3A3838','#808080','#DCDCDC'
)
names(chromHMM_palette) =chrom_states

##-----------
## functions
##-----------
## count SNPs
count_snps = function(x){x=x[,c(..range_keys)]%>%unique()%>%nrow()}

## Calculate MIAF/DAF
calculate_maf = function(x){
    df=copy(x)[
        ,MAF := ifelse(ancestry=='archaic', ## MIAF for aSNPs
            ifelse(POP_ARCH_REF > POP_ARCH_ALT,
            POP_ARCH_REF / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
            POP_ARCH_ALT / (POP_ARCH_REF + POP_ARCH_ALT + POP_NOTARCH_REF + POP_NOTARCH_ALT)
            ), ## DAF for naSNPs
            ifelse(ANC == 0, 
            POP_NOTARCH_ALT / (POP_ARCH_REF+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT),
            POP_NOTARCH_REF / (`POP_ARCH_REF`+POP_ARCH_ALT+POP_NOTARCH_REF+POP_NOTARCH_ALT))
            )
            ][
                ,state_allele := ifelse(ancestry=='archaic', ## label the state of the main introgressed archaic allele
                ifelse(
                    (ANC==0 & POP_ARCH_REF>POP_ARCH_ALT)|(ANC==1 & POP_ARCH_ALT>POP_ARCH_REF),'ancestral',
                    ifelse(POP_ARCH_REF==POP_ARCH_ALT,'not_determined','derived')
                    ),
                    'derived' ## all naSNPs that will be considered are derived. SNPs fixed for the ancestral allele are removed
                    )
                    ]
    return(df)

}


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
# calculate_odds_ratio = function(asnps,nasnps,merging_keys,fulljoin,numbsnps_col){
#     table = merge(asnps,nasnps,by=c(merging_keys),all=fulljoin)
#     table[is.na(table)] = 0
#     table = split(table,by=merging_keys)%>%
#     lapply(function(x)x=x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%
#     fisher.test()
#     )
#     numb_asnps <- copy(asnps)%>%setnames(old=c(numbsnps_col),new=c('numbsnps_element'))
#     numb_asnps <- numb_asnps[,c(..merging_keys,'numbsnps_element')]%>%split(by=merging_keys)%>%lapply(function(x)x[,c(merging_keys):=NULL])
      
#     fisher_test_results = copy(table)%>%
#     lapply(function(x)x=data.table(
#         'p'=x$p.value,
#         'odds_ratio'=x$estimate,
#         'lower_ci'=x$conf.int[[1]],
#         'upper_ci'=x$conf.int[[2]]
#     ))
#     fisher_test_results = Map(mutate,fisher_test_results,elements=names(fisher_test_results))%>%rbindlist()
#     numb_asnps = Map(mutate,numb_asnps,elements=names(numb_asnps))%>%rbindlist()

#     fisher_test_results = fisher_test_results[numb_asnps,on='elements',nomatch=0]

#     return(fisher_test_results)
# }

calculate_or <- function(snps_oi,snps_noi,merging_keys){
  oi_vs_noi_or <- merge(snps_oi,snps_noi,by=c(merging_keys))%>%
  dplyr::select(c(contains('numb'),contains(all_of(merging_keys))))

  fisher_test <- copy(oi_vs_noi_or)%>%split(by=c(merging_keys))%>%
  lapply(
    function(x){
    x <- x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%fisher.test()
    x <- data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
        )
    }
  )
  or_results <- Map(mutate,fisher_test,elements=names(fisher_test))%>%rbindlist()
  return(or_results)
}

## adjust pvalues
adjust_pvalues=function(x){
  df=copy(x)
  pvals_df=copy(df)
  pvals_df=pvals_df$p
  pvals_df_adjusted=p.adjust(pvals_df,'fdr')%>%as.data.table() %>% setnames('p_adj')
  pvals_df_adjusted=pvals_df_adjusted[
    ,p_signif:= ifelse(`p_adj`<=0.0001,'****',
                             ifelse(`p_adj`>0.0001 &`p_adj`<=0.001,'***',
                                    ifelse(`p_adj`>0.001 & `p_adj`<=0.01,'**',
                                           ifelse(`p_adj`>0.01 & `p_adj`<=0.05,'*',' '))))
    ]
   df_final=cbind(df,pvals_df_adjusted)
  return(df_final)
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
read_cres = function(dir,cells){
    cres = list.files(dir,recursive = F,full.names = T,patter='new_*') %>%
    lapply(
        function(y) z=fread(y,sep='\t',header = T)[
            chrom_state %in% c('1_TssA','2_TssAFlnk','3_TxFlnk',"6_EnhG","7_Enh")
            ][
            ,MAF:=round(MAF,2)
            ][
                ,chrom_state:=NULL
                ]%>%unique()
                # %>%setnames(old=12,new = 'allele_of_interest')
        )
    if(missing(cells)) {
        cres_allcells = copy(cres)
    return(cres_allcells)
    } else {
        cres_selected_cells = lapply(cres,function(x)x=x[cell_line %in% cells])
    return(cres_selected_cells)
    }
}

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

## create supp table
make_supp_table <- function(table,element_name){
    supp_table <- copy(table)%>%dplyr::select(
        c(contains('element'), contains('odds_ratio'),'lower_ci','upper_ci',contains('p'),'ancestry')
        )%>%setnames(old=c('elements','numbsnps_element'),new=c(element_name,'number_aSNPs'))
    return(supp_table)
}









