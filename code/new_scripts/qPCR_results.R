## use this script to analyse qPCR results
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)

options(width = 150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

qpcr_dir='./qPCR_results/'
table_dir='./Results/Tables/'
plot_dir='./Results/Plots/qPCR/'

##----------------------------------------
## first determine primer efficiencies 
##----------------------------------------
sample_code = fread(paste(qpcr_dir,'sample_names.txt',sep=''),sep='\t',header=T)

random_primers = c('11R','12R')
primer_eff_results=fread(paste(qpcr_dir,'primer_efficiencies/2021_07_05_summary_primer_efficiencies.txt',sep=''),sep='\t',header=T)
primer_eff_results = primer_eff_results[
    !Sample_Name %in% c('Stock','NTC')
    ][
        !Target_Name %in% random_primers
    ][
        ,Ct:=as.numeric(CT)
        ][
            ,CT:=NULL
]%>%na.omit()

primer_eff_results = primer_eff_results[
    ,log10_dilution:=ifelse(
        Sample_Name %in% '1/10',log10(1/10),
        ifelse(Sample_Name %in% '1/100',log10(1/100),
        ifelse(Sample_Name %in% '1/1000',log10(1/1000),log10(1/10000)
        )
        )
    )
    ][
        ,c('Sample_Name','Well_Position'):=NULL
]

## get mean Ct
primer_eff_mean_Ct = copy(primer_eff_results)%>%split(by=c('log10_dilution','Target_Name'))%>%lapply(function(x)x=x[,meanCt:=mean(Ct)][,Ct:=NULL]%>%unique())
primer_eff_mean_Ct = rbindlist(primer_eff_mean_Ct)%>%split(by='Target_Name')

## calculate primer efficiency
calculate_primer_efficiency = function(table){
    linear_model=lm(meanCt~log10_dilution,table)
    slope=linear_model$coefficients[[2]]
    efficiency=round(10^(-1/slope),2)
    return(efficiency)
}

primer_efficiency <- lapply(primer_eff_mean_Ct,function(x)calculate_primer_efficiency(x))%>%lapply(function(x)data.table(efficiency=x))
primer_efficiency <- Map(mutate,primer_efficiency,primer=names(primer_efficiency))%>%rbindlist()
primer_efficiency <- inner_join(primer_efficiency,sample_code,by=c('primer'='Code'))
primer_efficiency <- copy(primer_efficiency)[
    ,oligo:=paste(Oligo,allele_of_interest,sep='_')
]%>%setorderv('Oligo',1)

primer_colors <- c(
    "#001219","#005f73","#0a9396",
    "#94d2bd","#e9d8a6","#ee9b00",
    "#ca6702","#bb3e03","#ae2012",
    "#9b2226","#cb997e","#ddbea9",
    "#ffe8d6","#b7b7a4","#a5a58d",
    "#6b705c",'#f28482','#f6bd60'
)
names(primer_colors) = primer_efficiency$Oligo

## plot primer efficiencies
pdf(paste(plot_dir,'primer_efficiency_plot.pdf',sep=''),width=7,height=7)
ggplot(primer_efficiency,aes(x=oligo,y=efficiency,fill=Oligo))+
geom_bar(col='black',stat='identity')+
scale_fill_manual(values=primer_colors)+
geom_hline(yintercept=2,linetype='dashed',size=1)+
xlab('')+ylab('Efficiency')+
theme_bw()+
theme(
    axis.text.x=element_blank(),
    axis.ticks.x =element_blank()
    )
dev.off()

##-----------------
## now read results
##-----------------
bad_primers = c('1R','2R','9R','10R')
remove_snps = copy(primer_efficiency)[primer %in% bad_primers][,rsid:=sapply(strsplit(Oligo, "_\\s*"), tail, 1)]

read_qPCR_results=function(x){
    qPCR_results = fread(paste(qpcr_dir,x,sep=''),sep='\t',header=T)[
        !Sample_Name %in% 'NTC'
        ][
            ,c("Sample_Name",'Target_Name','CT')
            ][
                ,Ct:=as.numeric(CT)
                ][
                    ,CT:=NULL
                    ][
                        ,sample:=ifelse(Sample_Name %like% 'pDNA','plasmid','mRNA')
        ]%>%inner_join(sample_code,by=c('Target_Name'='Code'))

    qPCR_results = copy(qPCR_results)[
        ,replicate:=ifelse(Sample_Name %in% 'NTC','rep1',gsub('.*_','',Sample_Name))
        ][
            ,Sample_Name:=NULL
    ]
    return(qPCR_results)

}

png21_results = read_qPCR_results('/test_runs/2021_07_17/png21_qPCR_results.txt')[,cell_line:='png21']
png15_results = read_qPCR_results('/test_runs/2021_07_17/png15_qPCR_results.txt')[,cell_line:='png15']

combined_results =rbind(png21_results,png15_results)

combined_results <- combined_results[!Oligo %like% 'random']

## look at distribution differences between max and min Ct values between tech rep
qPCR_QC=copy(combined_results)%>%split(by=c('sample','Oligo','replicate','cell_line'))%>%lapply(
    function(x)x=x[
        ,Ct_diff:=max(Ct)-min(Ct)
        ][
            ,c('Ct_diff','sample','Oligo','cell_line')]%>%unique()
)%>%rbindlist()%>%setorderv('Oligo',1)

# pdf(paste(plot_dir,'QC_differences_Ct_values.pdf',sep=''),width=7,height=7)
# ggplot(qPCR_QC,aes(x=round(Ct_diff,1),fill=sample))+
# geom_bar(width = 0.05,position='dodge')+
# geom_vline(xintercept=0.5,linetype='dashed')+
# # facet_wrap(cell_line~.,ncol=1)+
# xlab('max(Ct)-min(Ct)')+
# theme(
#     panel.background =element_rect(fill = 'white', colour = 'black',size=1),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     strip.text.y = element_text(hjust = 0.5),
#     strip.background = element_rect(color = 'black', linetype = 'solid'),
#     strip.background.y = element_blank(),
#     strip.background.x =element_blank(),
#     legend.position='bottom',
#     legend.key = element_rect(fill = "white", colour = "black"),
#     axis.line = element_blank()
# )
# dev.off()

pdf(paste(plot_dir,'QC_differences_Ct_values_by_oligos.pdf',sep=''),width=7,height=7)
ggplot(qPCR_QC,aes(x=Oligo,y=Ct_diff,fill=Oligo))+
geom_boxplot()+
geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=0.05)+
scale_fill_manual(values=primer_colors)+
geom_hline(yintercept=0.5,linetype='dashed')+
facet_wrap(cell_line~.,ncol=1)+xlab(' ')+ylab('Ct difference')+
theme(
    panel.background =element_rect(fill = 'white', colour = 'black',size=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text.y = element_text(hjust = 0.5),
    strip.background = element_rect(color = 'black', linetype = 'solid'),
    strip.background.y = element_blank(),
    strip.background.x =element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank(),
    # legend.position='bottom',
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
)
dev.off()

## following this https://bitesizebio.com/24894/4-easy-steps-to-analyze-your-qpcr-data-using-double-delta-ct-analysis/

## GOI being tested = cDNA for allele of interest (AOI)
## housekeeper test = pDNA for AOI
## GOI control =  cDNA for other allele 
## housekeeper control = pDNA for other allele 

## dCt= GOI test/control - housekeeper test/control
## ddCt= dCt test - dCt control

## calculate mean Ct value per replicate
mean_Ct = copy(combined_results)%>%split(by=c('sample','Oligo','replicate','cell_line'))%>%lapply(
    function(x)x=x[
        ,meanCt:=mean(Ct)
        ][
            ,Ct:=NULL
            ][
                ,condition:=ifelse(Oligo %like% 'random','control','test')
                ]%>%unique()
)%>%rbindlist()
                    
mh_alleles = copy(mean_Ct)[allele_of_interest=='MH']
deni_alleles = copy(mean_Ct)[allele_of_interest=='Deni']
ptc_ref=copy(mean_Ct)[allele_of_interest=='ref']
ptc_alt=copy(mean_Ct)[allele_of_interest=='alt']
ntc = copy(mean_Ct)[Oligo=='random_15']

calculate_dCt =function(goi,housekeeper){
    goi=copy(goi)%>%inner_join(primer_efficiency,by=c('Target_Name'='primer','Oligo','allele_of_interest'))
    goi = copy(goi)[,n:=1:nrow(goi)] ## this is to make joining work otherwise it complains

    housekeeper=copy(housekeeper)%>%inner_join(primer_efficiency,by=c('Target_Name'='primer','Oligo','allele_of_interest'))
    housekeeper = copy(housekeeper)[,n:=1:nrow(housekeeper)] ## and this as well

    # deltaCt=purrr::map2(goi,housekeeper,function(x,y)x[
    #     y,on=c('replicate','cell_line'),nomatch=0
    #     ][
    #         ,dCt:=meanCt-i.meanCt
    #         ][
    #             ,E:=(efficiency^(-dCt))
    #             ][
    #                 ,c('Target_Name','Oligo','allele_of_interest','replicate','dCt','E','cell_line')
    #                 ]
    #     )%>%rbindlist()
    
    deltaCt=goi[
        housekeeper,on=c('replicate','cell_line','n'),nomatch=0
        ][
            ,dCt:=meanCt-i.meanCt
            ][
                ,E:=(efficiency^(-dCt))
                ][
                    ,c('Target_Name','Oligo','allele_of_interest','replicate','dCt','E','cell_line','n')
                    ]
    return(deltaCt)
}

calculate_ddCt =function(test,control){

    test_cDNA=copy(test)[sample=='mRNA']
    test_pDNA=copy(test)[sample!='mRNA']

    control_cDNA=copy(control)[sample=='mRNA']
    control_pDNA=copy(control)[sample!='mRNA']

    test_dCt=calculate_dCt(test_cDNA,test_pDNA)
    # test_dCt = test_dCt[,n:=1:nrow(test_dCt)] ## this is to make joining work otherwise it complains
    
    control_dCt=calculate_dCt(control_cDNA,control_pDNA)
    # control_dCt = control_dCt[,n:=1:nrow(control_dCt)] ## and this as well

    ddCt=test_dCt[
        control_dCt,on=c('n','replicate','cell_line'),nomatch=0
        ][
            ,ddCt:=dCt-i.dCt
            ][
                ,R:=E/i.E
            ]


    ddCt_and_fc = copy(ddCt)[,c('Oligo','allele_of_interest','R','ddCt','cell_line','replicate')][,rsid:=sub('.*\\_', '', Oligo)]%>%unique()

    return(ddCt_and_fc)
}

ptc_alt_ddCt=calculate_ddCt(ptc_alt,ptc_ref)
deni_ddCt=calculate_ddCt(deni_alleles,mh_alleles)

# nct_ddCt = calculate_ddCt(ntc[Oligo=='random_15'],ntc[Oligo=='random_22'])

# all_ddCts=list(ptc_alt_ddCt,deni_ddCt)
# all_ddCts =Map(mutate,all_ddCts,type=c('ptc','denisova'))%>%rbindlist()

all_ddCts=rbind(ptc_alt_ddCt,deni_ddCt)[!rsid %in% c(remove_snps$rsid)]

## perform stat test to see if relative expression is significantly different from 1
stat_test_expr=copy(all_ddCts)%>%split(by=c('rsid','cell_line'))%>%lapply(
    function(x)x=x[
        ,pval:=t.test(x$R,mu=1)$p.value
        ][
            ,p_signif:= ifelse(`pval`<=0.0001,'****',
                             ifelse(`pval`>0.0001 &`pval`<=0.001,'***',
                                    ifelse(`pval`>0.001 & `pval`<=0.01,'**',
                                           ifelse(`pval`>0.01 & `pval`<=0.05,'*',' '))))
    ][
        ,c('rsid','cell_line','pval','p_signif')
    ]%>%unique()
)%>%rbindlist()
# adjust_pvalues()

all_ddCts=all_ddCts[,c('Oligo','ddCt','allele_of_interest','replicate'):=NULL][, sapply(.SD, function(x) list(mean=mean(x), sd=sd(x))),by=.(cell_line,rsid)]
colnames(all_ddCts)[c(3:4)]=c('mean','sd')

all_ddCts=all_ddCts[stat_test_expr,on=c('rsid','cell_line'),nomatch=0]

rsid_cols <- c('#772F1A','#585123','#F2A65A','#EEC170','#F58549','#606c38','#EDE7E3')

pdf(paste(plot_dir,'oligo_fc.pdf',sep=''),width=4,height=5)
ggplot(all_ddCts, aes(x=rsid, y=mean,color=rsid,label = p_signif)) + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0, position=position_dodge(0.1))+
    # scale_color_manual(values=c(rep(c("#82c0cc"),5),"#ede7e3"),labels=c(rep('Denisova',5),'Tewhey et al 2016 \n positive control'))+
    scale_color_manual(values=c(rsid_cols),labels=c(rep('Denisova',5),'Tewhey et al 2016 \n positive control'))+
    geom_point(position=position_dodge(0.1))+
    xlab(' ')+ylab('Relative expression')+
    geom_hline(yintercept=1,linetype='dashed')+
    facet_wrap(cell_line~.,ncol=1)+
    geom_text(aes(y= all_ddCts$mean+all_ddCts$sd+0.001),size=5,color='black')+
    theme(
    panel.background =element_rect(fill = 'white', colour = 'black',size=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text.y = element_text(hjust = 0.5),
    strip.background = element_rect(color = 'black', linetype = 'solid'),
    strip.background.y = element_blank(),
    strip.background.x =element_blank(),
    legend.position='bottom',
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
    )
dev.off()

## check if relative expression correlates with delta PWM, i.e. higher affinities = higher expression 
## and/or allele frequency in PNG

## read deni files with all info for tfbs snps
denisova_tfbs_snps=fread(paste(table_dir,'supp_file_tfbs_snps_w_targets_Denisova.txt',sep=''),sep='\t',header=T)

deni_ddCt =deni_ddCt[
    denisova_tfbs_snps,on='rsid',nomatch=0
    ][
        ,c(..range_keys,'gene','R','cell_line','replicate','MAF','i.allele_of_interest','geneSymbol','alleleDiff','rsid')
        ]%>%split(by=c(range_keys))%>%lapply(
            function(x)x=x[
                ,avg_alleleDiff:=mean(alleleDiff)
                ]
)%>%rbindlist()

# ## plot distribution DPWMs for those SNPs
distribution_dPWM=copy(deni_ddCt)[,c('rsid','alleleDiff','geneSymbol')]%>%unique()

pdf(paste(plot_dir,'rsid_dPWM.pdf',sep=''),width=7,height=7)
ggplot(distribution_dPWM,aes(x=rsid,y=alleleDiff,fill=rsid))+
geom_boxplot()+
xlab(' ')+ylab('DPWM')+
geom_hline(yintercept=0,linetype='dashed')+
theme(
    panel.background =element_rect(fill = 'white', colour = 'black',size=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text.y = element_text(hjust = 0.5),
    strip.background = element_rect(color = 'black', linetype = 'solid'),
    strip.background.y = element_blank(),
    strip.background.x =element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
    )
dev.off()

