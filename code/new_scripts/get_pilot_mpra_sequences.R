library(data.table)
library(dplyr)
library(magrittr)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)


setwd('/data/projects/punim0586/dvespasiani/')
options(width = 150)


table_dir='./Archaic_introgression_in_PNG/Results/Tables/'


adapter_5 = 'ACTGGCCGCTTGACG'
adapter_3 = 'TAATCTCGGTCTGACTATCG'

downstream_seq_nobarcode = 'CACTGCGGCTCCTGCAGAGTACCTAGACGTCCATCGTACGCTTATCCGGACTGATGAGT'

barcodes = list(
  'gatcaggtgacactagcagg','ggacgaccttgttcgtgaac','cctaatgtcgacgcctattg',
  'caaactcgtattgggcctag','ttgttggtccccggaaccat','acgttcccggttagtcgaat',
  'ggcagatgtgaaagaaccgg','aggggtgacccctacaacat','ctggacggagtacggtaagg',
  'taggtgtgaaaccacttagc','cggacggtgttgtcctcata','tcgactctatagaggcatc',
  'ccatgtcgttactgaacgta','tcagacagatcacgtcatcg','atagagattacatctctaac',
  'ggtgatccgtcacgtatgag','aaggtatctgctcgacaatc','gctagactcagctgctcgac',
  'acgcccaggcatagttttag','gagtactagtcagtggcact'
)

## function to extract from fasta file the sequence of interest and remove the adapters
get_sequence_of_interest = function(x,soi){
    sequence = readDNAStringSet(x)%>%as.data.frame()%>%mutate('seqID'=rownames(.))%>%as.data.table()
    sequence = sequence[
        ,no_adapters:=substring(x, 16, 185)
        ][
            seqID %like% soi
            ][
                ,x:=NULL
                ]
    return(sequence)

}

get_oligo = function(x){
    df=copy(x)
    df=df[
        ,oligo:=as.character(getSeq(Hsapiens,df$seqnames,start = df$start-84,end =df$start+85))
    ]
    substr(df$oligo, 85,85)=df$allele

    return(df)
}


random_sequences = get_sequence_of_interest('./Generate_random_DNA_sequences/output/final_random_sequences.fa','random_15$|random_22$')
tewhey_control_sequences = get_sequence_of_interest('./Files/Functional_validation_denisova_snps/Sequences_MPRA/full_library/fasta_files/controls/control_snp_sequences.fa','rs9283753')

tewhey_control_sequences = tewhey_control_sequences[
    ,c("seqnames", "start",'rsid','allele') := tstrsplit(seqID, "_", fixed=TRUE)
    ][
        ,c('seqID','no_adapters'):=NULL
        ][
            ,start:=as.numeric(start)
]%>%get_oligo()
tewhey_control_sequences = tewhey_control_sequences[,seqID:=paste(seqnames,start,allele,rsid,sep='_')]

oas_snps = fread('./Archaic_introgression_in_PNG/Results/Tables/OAS_snps.txt',sep='\t',header=T)
oas_snps_ref = copy(oas_snps)[,c(1:4,11)][,state_allele:='ref']%>%setnames(old='REF',new='allele')
oas_snps_alt = copy(oas_snps)[,c(1:3,5,11)][,state_allele:='alt']%>%setnames(old='ALT',new='allele')

oas_snps_oligos = rbind(oas_snps_ref,oas_snps_alt)%>%get_oligo()
oas_snps_oligos = oas_snps_oligos[,seqID:=paste(seqnames,start,allele,rsid,sep='_')]

all_oligos = rbind(
    oas_snps_oligos[,c('seqID','oligo')],
    tewhey_control_sequences[,c('seqID','oligo')],
    copy(random_sequences)%>%setnames(old='no_adapters',new='oligo')
)

## add adapters, barcodes and stuff
all_oligos = all_oligos[
    ,final_oligo:=paste(adapter_5,oligo,downstream_seq_nobarcode,sep='')
]%>%split(by='seqID')

all_oligos = purrr::map2(all_oligos,barcodes,function(x,y)x=x[,final_oligo:=paste(final_oligo,y,adapter_3,sep='')])%>%rbindlist()


## check for enzymes
enzyme_free=function(x,column_to_check){
  BsiWI='CGTACG'
  SfiI='GGCC[ATCG]{1,5}GGCC'
  

  match_enzyme=function(y,enzyme){
    df=copy(y)
    name_enzy = deparse(substitute(enzyme))
    df=df[
      ,match_enzyme:=ifelse(grepl(enzyme, dna_seq, ignore.case=TRUE),enzyme,'n')
      ][
          ,which_enzyme:= name_enzy
          ][
              ,number_matches:= stringi::stri_count_regex(df$dna_seq, enzyme)
          ]
    return(df)
  }
  
  df_check = copy(x)[,dna_seq:=column_to_check]
  BsiWI_results = match_enzyme(df_check,BsiWI)
  SfiI_results = match_enzyme(df_check,SfiI) 

  enzyme_results = rbind(BsiWI_results,SfiI_results)[,dna_seq:=NULL]

  return(enzyme_results)
}

all_oligos_enyzmefree = enzyme_free(all_oligos,all_oligos$final_oligo)
all_oligos_enyzmefree = all_oligos_enyzmefree[match_enzyme=='n'][,c('seqID','final_oligo')]
colnames(all_oligos_enyzmefree)= c('Name',"Sequence(5-3)")


## get primers
R_primers_barcodes = copy(barcodes)%>%lapply(function(x)DNAString(x)%>%reverseComplement()%>%as.character())
names(R_primers_barcodes) = all_oligos_enyzmefree$Name

F_primer_1st_pcr = adapter_5
R_primer_1st_pcr = DNAString(adapter_3)%>%reverseComplement()%>%as.character()

F_primer_2nd_pcr = paste('GCTCAGAACAGCACATCTGGCCTA',adapter_5,sep='')
R_primer_2nd_pcr = paste('GTTTAAGGCCTCCGTGGCC',R_primer_1st_pcr,sep='')

pcr_primers = list(F_primer_1st_pcr,R_primer_1st_pcr,F_primer_2nd_pcr,R_primer_2nd_pcr)
names(pcr_primers) =c('F_primer_1st_pcr','R_primer_1st_pcr','F_primer_2nd_pcr','R_primer_2nd_pcr')

all_primers = c(pcr_primers,R_primers_barcodes)%>%lapply(function(x)data.table('Sequence(5-3)'=x))
all_primers = Map(mutate,all_primers,Name=names(all_primers))%>%rbindlist()

sequences_pilot_mpra = list(all_oligos_enyzmefree,all_primers)
names(sequences_pilot_mpra) = c('oligos','primers')

## write file
write.xlsx(sequences_pilot_mpra,paste(table_dir,'oligo_sequences.xlsx',sep=''),append=T)



