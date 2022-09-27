## use this script to plot regions with Gviz
## ADD DHS info
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(rtracklayer)
library(Gviz)
library(biomaRt)
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocParallel)


options(width=150)
setwd('/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

chrom_state_dir = '../Annotation_and_other_files/Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states'
plot_dir = './Results/Plots/GenomeBrowser/'
table_dir='./Results/Tables/'

human_mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")

##---------------
## read files
##---------------
## OAS snps
denisova_oas_snps=fread(paste(table_dir,'supp_file_tfbs_snps_w_targets_Denisova.txt',sep=''),sep='\t',header=T)[gene%like% 'OAS']

pwd_with_strongest_effect=copy(denisova_oas_snps)[, .SD[which.max(abs(alleleDiff))], by=.(seqnames,start,end,REF,ALT)]$providerName

# ## rerun motifbreak for these 8 SNPs because u need to get the PWM
# hocomoco_cores = c(
#   'HOCOMOCOv11-core-A','HOCOMOCOv11-core-B','HOCOMOCOv11-core-C',
#   'HOCOMOCOv11-secondary-A','HOCOMOCOv11-secondary-B','HOCOMOCOv11-secondary-C'
# )

# jaspar2018 = subset(MotifDb, organism=='Hsapiens' & dataSource=='jaspar2018')
# hocomoco = subset(MotifDb, organism=='Hsapiens' & dataSource %in% hocomoco_cores)

# ## read in oas results in motifbreak-ready format
# oas_snps_motifbreak_input= snps.from.file(
#     file = './Motifbreak/input_files/denisova.bed',
#     search.genome = BSgenome.Hsapiens.UCSC.hg19,
#     format = "bed"
# )

# oas_snps_motifbreak_input=oas_snps_motifbreak_input[end(oas_snps_motifbreak_input)%in% denisova_oas_snps$start]

# oas_snps_motifbreak_output=motifbreakR(
#     snpList = oas_snps_motifbreak_input, filterp = TRUE,
#     pwmList = c(hocomoco,jaspar2018),threshold = 1e-5,
#     method = "log",bkg = c(A=0.3, C=0.2, G=0.2, T=0.3), 
#     BPPARAM = BiocParallel::bpparam()
# )

# pwms_to_plot=copy(oas_snps_motifbreak_output)
# pwms_to_plot=pwms_to_plot[pwms_to_plot$providerName %in% pwd_with_strongest_effect]
# pwms_to_plot=split(pwms_to_plot,start(pwms_to_plot))

# names(pwms_to_plot)=oas_snps_motifbreak_input$SNP_id

# ## plot PWMs
# pdf(paste(plot_dir,'motifbreak_PWMs.pdf',sep=''),width=7,height=5)
# lapply(pwms_to_plot,function(x)plotMB(x,names(x)))
# dev.off()


## Roadmap chrom states immune cells
immune_cells_chromstates=list.files(chrom_state_dir,recursive=F,full.names=T,pattern="Bc|Tc")%>%lapply(function(x)fread(x,sep=' ',header=T))

##------------------------
## Create OAS gene track
##------------------------
oas_gene_region=copy(denisova_oas_snps)[,start:=min(start)-7000][,end:=max(end)+7000][,c(..range_keys)]%>%unique()

oas_gtrack = GenomeAxisTrack(
    chromosome = unique(as.character(oas_gene_region$seqnames)),
    start = min(oas_gene_region$start),
    end =  max(oas_gene_region$end),
    size = 0.5
)

biomart_track = BiomartGeneRegionTrack(
    chromosome = as.character(unique(oas_gene_region$seqnames)),
    start = as.numeric(oas_gene_region$start),
    end = as.numeric(oas_gene_region$end),
    biomart =human_mart,
    genome='hg19',
    name='Target gene',
    size=0.2,
    min.height = 0.5 
)

biomart_track=biomart_track[biomart_track@range$symbol %in% c('OAS2','OAS3')]

pdf(paste0(plot_dir,'oas_track.pdf',sep=''),width = 10, height = 7)
plotTracks(
    c(oas_gtrack,biomart_track),
    from = oas_gene_region$start,
    to =  oas_gene_region$end,
    collapseTranscripts='longest', 
    transcriptAnnotation = 'symbol'
)
dev.off()


## create 2 SNPs tracks, 1 per OAS target gene
get_snp_region=function(snp){
    region=copy(snp)[
            ,start:=min(start)-40000
            ][
                ,end:=max(end)+40000
                ][
                    ,c(..range_keys)
                    ]%>%unique()
    return(region)
}
get_chromstate = function(region,chromstate){
    grange = copy(region)
    chromstates = copy(chromstate)[
        seqnames %in% unique(grange$seqnames)
        ][
            start>=unique(grange$start)+1
            ][
                end<=unique(grange$end)
                ][
                    ,state_numb:=as.numeric(gsub('\\_.*','',chrom_state))
                    ][
                        ,c('chrom_state','cell_type'):=NULL
                        ]%>%unique()%>%makeGRangesFromDataFrame(keep.extra.columns=T)
    return(chromstates)

}
# get_chromstate = function(region,chromstate){
#     grange = copy(region)
#     setkeyv(grange,range_keys)
#     chromstates = copy(chromstate)
#     setkeyv(chromstates,range_keys)

#     palette=data.table(palette=chromHMM_palette,chrom_state=names(chromHMM_palette))
#     overlaps=foverlaps(grange,chromstates,type='any')[palette,on='chrom_state',nomatch=0][
#         ,c(..range_keys,'palette')]%>%unique()%>%
#         makeGRangesFromDataFrame(keep.extra.columns=T)

#     return(overlaps)

# }
get_snp_data=function(x){
    df=copy(x)
    snps=copy(x)[,c(..range_keys,'MAF')]%>%unique()
    region=get_snp_region(df)
    chromstate_bcells=get_chromstate(region,immune_cells_chromstates[[1]])
    chromstate_tcells=get_chromstate(region,immune_cells_chromstates[[2]])

    snp_pwm_scores=copy(df)[providerName%in%pwd_with_strongest_effect][,c(..range_keys,'alleleDiff')]

    snp_data=list(snps,region,chromstate_bcells,chromstate_tcells,snp_pwm_scores)
    names(snp_data)=c('snps','region','chromstate_bcells','chromstate_tcells','pwm_scores')
    return(snp_data)
}

oas2_snps_data=get_snp_data(denisova_oas_snps[gene=='OAS2'])
oas3_snps_data=get_snp_data(denisova_oas_snps[gene=='OAS3'])

## Roadmap DHS immune cells
immune_cells_bigwigs = list.files('../Annotation_and_other_files/Roadmap_data/immune_cells_DNase_bigwig',recursive=T,full.names=T,pattern= 'fc.signal.bigwig')
##-----------------
## Create tracks
##-----------------
get_bigwig_track = function(file,region,colors){
    track_colors = as.list(c(rep(colors,length(file))))
    
    files = copy(file)
    sample_names = lapply(files,function(z)gsub('\\_.*','',basename(z)))
    sample_names = lapply(sample_names,function(x)gsub('\\..*','',x))

   track = purrr::pmap(
        list(files,sample_names,colors),
        function(x,y,z)
        DataTrack(
            range = x, 
            genome = "hg19",
            type = "h", 
            chromosome = unique(as.character(region$seqnames)),
            start = unique(region$start),
            end = unique(region$end),
            name = y,
            col=z,
            col.histogram=z,
            size=0.2
            )
        )
    track=OverlayTrack(track)
    return(track)
}

get_genome_browser_tracks=function(data){

    oas_gtrack = GenomeAxisTrack(
        chromosome = unique(as.character(data$snps$seqnames)),
        start = min(data$snps$start),
        end =  max(data$snps$end),
        size = 0.5
        )
    names(chromHMM_palette)=NULL
    bcells_chromstate_track = DataTrack(
            data$chromstate_bcells,
            type = "heatmap",
            name = "B Cells chrom state",
            gradient = chromHMM_palette,
            size = 0.2,
            min.height = 0.5
    )
    tcells_chromstate_track = DataTrack(
            data$chromstate_tcells,
            type = "heatmap",
            name = "T Cells chrom state",
            gradient = chromHMM_palette,
            size = 0.2,
            min.height = 0.5
    )

    # bcells_chromstate_track = AnnotationTrack(
    #     range=makeGRangesFromDataFrame(data$chromstate_bcells,keep.extra.columns=T),
    #     genome = 'hg19',
    #     fill=chromHMM_palette,
    #     col='black',
    #     name = "bcells_chromstate"
    # )
  
    # tcells_chromstate_track = AnnotationTrack(
    #     range=makeGRangesFromDataFrame(data$chromstate_tcells,keep.extra.columns=T),
    #     genome = 'hg19',
    #     fill=chromHMM_palette,
    #     col='black',
    #     name = "tcells_chromstate"
    # )

    snpfreq_track = DataTrack(
            range = makeGRangesFromDataFrame(unique(data$snps),keep.extra.columns=T),
            name = "MIAF",
            type='p',
            size = 0.3
    )

    pwm_track = DataTrack(
            range = makeGRangesFromDataFrame(data$pwm_scores,keep.extra.columns=T),
            name = "DPWM",
            type='h',
            size = 0.3,
            lwd = 0.3,
            min.height = 0.3
    )
    
    bcells_dhs_track = get_bigwig_track(immune_cells_bigwigs[c(1:5)],data$region,'palegreen4')
    tcells_dhs_track = get_bigwig_track(immune_cells_bigwigs[c(6:7)],data$region,'chartreuse4')

    all_tracks = c(oas_gtrack,bcells_dhs_track,tcells_dhs_track,snpfreq_track,pwm_track,bcells_chromstate_track,tcells_chromstate_track)

    return(all_tracks)
}

## plot everything
oas2_all_tracks = get_genome_browser_tracks(oas2_snps_data)
oas3_all_tracks = get_genome_browser_tracks(oas3_snps_data)

pdf(paste0(plot_dir,'oas2_genome_browser.pdf',sep=''),width = 10, height = 7)
plotTracks(
    oas2_all_tracks,
    from=oas2_snps_data$region$start,
    to=oas2_snps_data$region$end
    # stacking = "dense"
)
dev.off()

pdf(paste0(plot_dir,'oas3_genome_browser.pdf',sep=''),width = 10, height = 7)
plotTracks(
    oas3_all_tracks,
    from=oas3_snps_data$region$start,
    to=oas3_snps_data$region$end
    # stacking = "dense"
)
dev.off()

# oas_gtrack = GenomeAxisTrack(
#     chromosome = unique(as.character(denisova_oas_snps$seqnames)),
#     start = min(denisova_oas_snps$start),
#     end =  max(denisova_oas_snps$end),
#     size = 0.5
# )

# ## now for snps
# snps_track = GenomeAxisTrack(
#     chromosome = unique(as.character(oas_gene_region$seqnames)),
#     start = min(oas_gene_region$start),
#     end =  max(oas_gene_region$end),
#     size = 0.5
# )


# names(chromHMM_palette)=NULL ## only for here

# bcells_chromstate_track = DataTrack(
#         oas_chromstate_bcells,
#         type = "heatmap",
#         name = "B Cells chrom state",
#         gradient = chromHMM_palette,
#         size = 0.2,
#         min.height = 0.5
# )

# tcells_chromstate_track = DataTrack(
#         oas_chromstate_tcells,
#         type = "heatmap",
#         name = "T Cells chrom state",
#         gradient = chromHMM_palette,
#         size = 0.2,
#         min.height = 0.5
# )


# snpfreq_track = DataTrack(
#         range = makeGRangesFromDataFrame(unique(denisova_oas_snps[,c(..range_keys,'MAF')]),keep.extra.columns=T),
#         name = "MIAF",
#         type='p',
#         size = 0.3
# )

# pwm_track = DataTrack(
#         range = makeGRangesFromDataFrame(snp_pwm_scores,keep.extra.columns=T),
#         name = "DPWM",
#         type='boxplot',
#         type='h',
#         size = 0.3,
#         lwd = 0.3,
#         min.height = 0.3
# )

# get_track = function(file,region,colors){
#     track_colors = as.list(c(rep(colors,length(file))))
    
#     files = copy(file)
#     sample_names = lapply(files,function(z)gsub('\\_.*','',basename(z)))
#     sample_names = lapply(sample_names,function(x)gsub('\\..*','',x))

#    track = purrr::pmap(
#         list(files,sample_names,colors),
#         function(x,y,z)
#         DataTrack(
#             range = x, 
#             genome = "hg19",
#             type = "h", 
#             chromosome = unique(as.character(region$seqnames)),
#             start = unique(region$start),
#             end = unique(region$end),
#             name = y,
#             col=z,
#             col.histogram=z,
#             size=0.2
#             )
#         )
#     track=OverlayTrack(track)
#     return(track)
# }
# bcells_track = get_track(immune_cells_bigwigs[c(1:5)],region_around_oas_snps,'palegreen4')
# tcells_track = get_track(immune_cells_bigwigs[c(6:7)],region_around_oas_snps,'chartreuse4')

## plot everything
all_tracks = c(gtrack,bcells_track,tcells_track,snpfreq_track,pwm_track,bcells_chromstate_track,tcells_chromstate_track,biomart_track)

pdf(paste0(plot_dir,'oas_genome_browser.pdf',sep=''),width = 10, height = 7)
plotTracks(
    all_tracks,
    from = region_around_oas_snps$start,
    to =  region_around_oas_snps$end,
    collapseTranscripts='longest', 
    transcriptAnnotation = 'symbol'
)
dev.off()