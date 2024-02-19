library(tidyverse)
library(IRanges)
library(dplyr)
library(plyranges)
library(ChIPpeakAnno)




#Obtain TSS coordinates for all human genome:

TSS <- TSS.human.GRCh38 %>% 
  addchr() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  dplyr::rename(Geneid=1) 

library(biomaRt)

Geneid <- TSS$Geneid

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

listAttributes(ensembl)


TSS_symbol <- getBM(attributes = c("external_gene_name","ensembl_gene_id"),filters = "ensembl_gene_id",values = Geneid,mart = ensembl) %>% 
  dplyr::rename(SYMBOL=1,
                Geneid=2)

TSS <- inner_join(TSS,TSS_symbol) %>% 
  dplyr::rename(chr=2,start=3,end=4) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

TSS  <- TSS[seqnames(TSS) %in% c("chr1","chr2","chr4","chr3","chr5","chr6","chr7",
                                 "chr8","chr9","chr10","chr11","chr12","chr13",
                                 "chr14","chr15","chr16","chr17","chr18","chr19",
                                 "chr20","chr21","chr22","chrX","chrY")] %>% 
  as.data.frame()%>% 
  dplyr::rename(chr_T=1,start_T=2,end_T=3, strand_T=5)

df <- data.frame(unique(TSS$seqnames))


#Calculate enhancers-promoters contact in  using all Neu D868D all peaks 


#Associating distal (>50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868D <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations (b6236370)/mega/NPC_Neuron-SETBP1_D868D/annotation.bed.gz", #Table of binned enhancers-promoters obtained with a costum script from a collaborator
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=1,
                start=2,
                end=3,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868D_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868D,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868D <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  dplyr::filter(Neu_D868D_peaks==1| Neu_D868D_up_dev==1) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868D_HiC <- join_overlap_inner(ATAC_peaks_neurons_D868D,Promoters_enhancers_interaction_Neu_D868D_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868D_HiC$x <- apply( ATAC_peaks_neurons_D868D_HiC[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868D_HiC  <- ATAC_peaks_neurons_D868D_HiC %>% 
  dplyr::rename(peaks=18)


ATAC_peaks_neurons_D868D_HiC <- ATAC_peaks_neurons_D868D_HiC[with(ATAC_peaks_neurons_D868D_HiC, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868D_HiC <- ATAC_peaks_neurons_D868D_HiC %>% 
  dplyr::rename(chr=1,
                chr2=9,
                start2=10,
                SYMBOL=8,
                end2=11)

#Associating proximal (<50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868D <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations/mega (2f06dcec)/NPC_Neuron-SETBP1_D868D/straw_norm.tsv.gz", #Table of binned enhancers-promoters obtained with a costum script from a collaborator
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=7,
                start=8,
                end=9,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868D_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868D,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868D <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  dplyr::filter(Neu_D868D_peaks==1| Neu_D868D_up_dev==1) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868D_HiC_p <- join_overlap_inner(ATAC_peaks_neurons_D868D,Promoters_enhancers_interaction_Neu_D868D_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868D_HiC_p$x <- apply(ATAC_peaks_neurons_D868D_HiC_p[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868D_HiC_p  <- ATAC_peaks_neurons_D868D_HiC_p %>% 
  dplyr::rename(peaks=18)


ATAC_peaks_neurons_D868D_HiC_p <- ATAC_peaks_neurons_D868D_HiC_p[with(ATAC_peaks_neurons_D868D_HiC_p, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868D_HiC_p <- ATAC_peaks_neurons_D868D_HiC_p %>% 
  dplyr::rename(chr=1,
                chr2=6,
                start2=7,
                SYMBOL=11,
                end2=8) 

ATAC_peaks_neurons_D868D_HiC <- rbind.data.frame(ATAC_peaks_neurons_D868D_HiC,ATAC_peaks_neurons_D868D_HiC_p) # Create a unique list of contact proximal and distal 

ATAC_peaks_neurons_D868D_HiC_plus <- ATAC_peaks_neurons_D868D_HiC %>%   # Modifying dataset to have the correct TSS site for genes on minus strand 
  dplyr::filter(strand_T==c("+")) %>% 
  dplyr::mutate(end_T=start_T+1000) %>% 
  dplyr::mutate(distanceToTSS=start_T-end) 

ATAC_peaks_neurons_D868D_HiC_minus <- ATAC_peaks_neurons_D868D_HiC %>% 
  dplyr::filter(strand_T==c("-")) %>% 
  dplyr::mutate(start_T=end_T-1000) %>% 
  dplyr::mutate(distanceToTSS=end_T-end) 



ATAC_peaks_neurons_D868D_HiC <- rbind.data.frame(ATAC_peaks_neurons_D868D_HiC_plus,ATAC_peaks_neurons_D868D_HiC_minus)

write_delim(ATAC_peaks_neurons_D868D_HiC,"SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_HiC",
            delim="\t", col_names = T)

ATAC_peaks_neurons_D868D_HiC.bedpe <- ATAC_peaks_neurons_D868D_HiC %>%    #Create tacks to visualize contact on IGV 
  dplyr::select(chr,start,end,chr_T,start_T,end_T,IF) %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_HiC.bedpe",
              delim="\t",col_names = F)


#Calculate enhancers-promoters contact using all Neu D868N all peaks


#Associating distal (>50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868N <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations (b6236370)/mega/NPC_Neuron-SETBP1_D868N/annotation.bed.gz", #Table of binned enhancers-promoters obtained with a costum script from a collaborator
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=1,
                start=2,
                end=3,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868N_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868N,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868N <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  dplyr::filter(Neu_D868N_peaks==1| Neu_D868N_up_dev==1) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868N_HiC <- join_overlap_inner(ATAC_peaks_neurons_D868N,Promoters_enhancers_interaction_Neu_D868N_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868N_HiC$x <- apply( ATAC_peaks_neurons_D868N_HiC[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868N_HiC  <- ATAC_peaks_neurons_D868N_HiC %>% 
  dplyr::rename(peaks=18)


ATAC_peaks_neurons_D868N_HiC <- ATAC_peaks_neurons_D868N_HiC[with(ATAC_peaks_neurons_D868N_HiC, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868N_HiC <- ATAC_peaks_neurons_D868N_HiC %>% 
  dplyr::rename(chr=1,
                chr2=9,
                start2=10,
                SYMBOL=8,
                end2=11)

#Associating proximal (<50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868N <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations/mega (2f06dcec)/NPC_Neuron-SETBP1_D868N/straw_norm.tsv.gz", #Table of binned enhancers-promoters obtained with a costum script from a collaborator
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=7,
                start=8,
                end=9,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868N_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868N,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868N <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  dplyr::filter(Neu_D868N_peaks==1| Neu_D868N_up_dev==1) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868N_HiC_p <- join_overlap_inner(ATAC_peaks_neurons_D868N,Promoters_enhancers_interaction_Neu_D868N_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868N_HiC_p$x <- apply(ATAC_peaks_neurons_D868N_HiC_p[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868N_HiC_p  <- ATAC_peaks_neurons_D868N_HiC_p %>% 
  dplyr::rename(peaks=18)


ATAC_peaks_neurons_D868N_HiC_p <- ATAC_peaks_neurons_D868N_HiC_p[with(ATAC_peaks_neurons_D868N_HiC_p, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868N_HiC_p <- ATAC_peaks_neurons_D868N_HiC_p %>% 
  dplyr::rename(chr=1,
                chr2=6,
                start2=7,
                SYMBOL=11,
                end2=8) 

ATAC_peaks_neurons_D868N_HiC <- rbind.data.frame(ATAC_peaks_neurons_D868N_HiC,ATAC_peaks_neurons_D868N_HiC_p) # Create a unique list of contact proximal and distal 

ATAC_peaks_neurons_D868N_HiC_plus <- ATAC_peaks_neurons_D868N_HiC %>%   # Modifying dataset to have the correct TSS site for genes on minus strand 
  dplyr::filter(strand_T==c("+")) %>% 
  dplyr::mutate(end_T=start_T+1000) %>% 
  dplyr::mutate(distanceToTSS=start_T-end) 

ATAC_peaks_neurons_D868N_HiC_minus <- ATAC_peaks_neurons_D868N_HiC %>% 
  dplyr::filter(strand_T==c("-")) %>% 
  dplyr::mutate(start_T=end_T-1000) %>% 
  dplyr::mutate(distanceToTSS=end_T-end) 



ATAC_peaks_neurons_D868N_HiC <- rbind.data.frame(ATAC_peaks_neurons_D868N_HiC_plus,ATAC_peaks_neurons_D868N_HiC_minus)

write_delim(ATAC_peaks_neurons_D868N_HiC,"SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_HiC",
            delim="\t", col_names = T)

ATAC_peaks_neurons_D868N_HiC.bedpe <- ATAC_peaks_neurons_D868N_HiC %>%    #Create tacks to visualize contact on IGV 
  dplyr::select(chr,start,end,chr_T,start_T,end_T,IF) %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_HiC.bedpe",
              delim="\t",col_names = F)

#Calculate enhancers-promoters contact using all Neu D868D peaks upregulated in Neuronal diff.

#Associating distal (>50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868D <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations (b6236370)/mega/NPC_Neuron-SETBP1_D868D/annotation.bed.gz",
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=1,
                start=2,
                end=3,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868D_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868D,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

ATAC_peaks_neurons_D868D_UP <- read_tsv("SETBP1_epigenomics/pipeline/Deseq2/Neu_D868DvsNPC_D868D_UP.tsv")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868D_HiC_up_dev <- join_overlap_inner(ATAC_peaks_neurons_D868D_UP,Promoters_enhancers_interaction_Neu_D868D_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868D_HiC_up_dev$x <- apply( ATAC_peaks_neurons_D868D_HiC_up_dev[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868D_HiC_up_dev <- ATAC_peaks_neurons_D868D_HiC_up_dev %>% 
  dplyr::rename(peaks=17)


ATAC_peaks_neurons_D868D_HiC_up_dev <- ATAC_peaks_neurons_D868D_HiC_up_dev[with(ATAC_peaks_neurons_D868D_HiC_up_dev, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868D_HiC_up_dev <- ATAC_peaks_neurons_D868D_HiC_up_dev %>% 
  dplyr::rename(chr=1,
                chr2=9,
                start2=10,
                SYMBOL=8,
                end2=11) %>% 
  dplyr::mutate(distanceToTSS=end-start_T)

#Associating proximal (<50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868D <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations (b6236370)/mega/NPC_Neuron-SETBP1_D868D/annotation.bed.gz",
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=7,
                start=8,
                end=9,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868D_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868D,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868D_HiC_up_dev_p <- join_overlap_inner(ATAC_peaks_neurons_D868D_UP,Promoters_enhancers_interaction_Neu_D868D_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868D_HiC_up_dev_p$x <- apply( ATAC_peaks_neurons_D868D_HiC_up_dev_p[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868D_HiC_up_dev_p <- ATAC_peaks_neurons_D868D_HiC_up_dev_p %>% 
  dplyr::rename(peaks=17)


ATAC_peaks_neurons_D868D_HiC_up_dev_p <- ATAC_peaks_neurons_D868D_HiC_up_dev_p[with(ATAC_peaks_neurons_D868D_HiC_up_dev_p, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868D_HiC_up_dev_p <- ATAC_peaks_neurons_D868D_HiC_up_dev_p %>% 
  dplyr::rename(chr=1,
                chr2=6,
                start2=7,
                SYMBOL=11,
                end2=8) %>% 
  dplyr::mutate(distanceToTSS=end-start_T)

ATAC_peaks_neurons_D868D_HiC_up_dev <- rbind.data.frame(ATAC_peaks_neurons_D868D_HiC_up_dev,ATAC_peaks_neurons_D868D_HiC_up_dev_p)
ATAC_peaks_neurons_D868D_HiC_up_dev <- ATAC_peaks_neurons_D868D_HiC_up_dev[with(ATAC_peaks_neurons_D868D_HiC_up_dev, ave(IF, peaks, FUN= max)==IF),]


write_delim(ATAC_peaks_neurons_D868D_HiC_up_dev,"SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_UP_HiC_dev",
            delim="\t", col_names = T)

ATAC_peaks_neurons_D868D_HiC_up_dev.bedpe <- ATAC_peaks_neurons_D868D_HiC_up_dev %>% 
  dplyr::select(chr,start,end,chr2,start2,end2,IF) %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_UP_HiC_dev.bedpe",
              delim="\t",col_names = F)

ATAC_peaks_neurons_D868D_HiC_up_dev_symbol <- ATAC_peaks_neurons_D868D_HiC_up_dev %>% # Extracting associated gene list
  dplyr::rename(SYMBOL=8) %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_UP_HiC_dev_SYMBOL",
              delim="\t",col_names = T)

#Calculate enhancers-promoters contact using all Neu D868N peaks upregulated in Neuronal diff.

#Associating distal (>50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868N <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations/mega/NPC_Neuron-SETBP1_D868N/annotation.bed.gz",
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=1,
                start=2,
                end=3,
                IF=5,
                SYMBOL=6)

Promoters_enhancers_interaction_Neu_D868N_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868N,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

ATAC_peaks_neurons_D868N_UP <- read_tsv("SETBP1_epigenomics/pipeline/Deseq2/Neu_D868NvsNPC_D868N_UP.tsv")%>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame()

ATAC_peaks_neurons_D868N_HiC_up_dev <- join_overlap_inner(ATAC_peaks_neurons_D868N_UP,Promoters_enhancers_interaction_Neu_D868N_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868N_HiC_up_dev$x <- apply( ATAC_peaks_neurons_D868N_HiC_up_dev[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868N_HiC_up_dev <- ATAC_peaks_neurons_D868N_HiC_up_dev %>% 
  dplyr::rename(peaks=17)


ATAC_peaks_neurons_D868N_HiC_up_dev <- ATAC_peaks_neurons_D868N_HiC_up_dev[with(ATAC_peaks_neurons_D868N_HiC_up_dev, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868N_HiC_up_dev <- ATAC_peaks_neurons_D868N_HiC_up_dev %>% 
  dplyr::rename(chr=1,
                chr2=9,
                start2=10,
                SYMBOL=8,
                end2=11) %>% 
  dplyr::mutate(distanceToTSS=start-start2)

#Associating proximal (<50 kb) peaks to a cognate promoters 

Promoters_enhancers_interaction_Neu_D868N <- read_delim("Setbp1_Gdrive/setbp1/pipeline/annotations (b6236370)/mega/NPC_Neuron-SETBP1_D868N/annotation.bed.gz",
                                                        delim="\t",col_names = F) %>%
  dplyr::rename(chr=7,
                start=8,
                end=9,
                IF=5,
                SYMBOL=6) 

Promoters_enhancers_interaction_Neu_D868N_TSS <- inner_join(Promoters_enhancers_interaction_Neu_D868N,TSS)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


ATAC_peaks_neurons_D868N_HiC_up_dev_p <- join_overlap_inner(ATAC_peaks_neurons_D868N_UP,Promoters_enhancers_interaction_Neu_D868N_TSS) %>% 
  as.data.frame()



paste_noNA <- function(x,sep=",") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }

sep="_"
ATAC_peaks_neurons_D868N_HiC_up_dev_p$x <- apply( ATAC_peaks_neurons_D868N_HiC_up_dev_p[ , c(1:3) ] , 1 , paste_noNA , sep=sep)
df

ATAC_peaks_neurons_D868N_HiC_up_dev_p <- ATAC_peaks_neurons_D868N_HiC_up_dev_p %>% 
  dplyr::rename(peaks=17)


ATAC_peaks_neurons_D868N_HiC_up_dev_p <- ATAC_peaks_neurons_D868N_HiC_up_dev_p[with(ATAC_peaks_neurons_D868N_HiC_up_dev_p, ave(IF, peaks, FUN= max)==IF),]

ATAC_peaks_neurons_D868N_HiC_up_dev_p <- ATAC_peaks_neurons_D868N_HiC_up_dev_p %>% 
  dplyr::rename(chr=1,
                chr2=6,
                start2=7,
                SYMBOL=11,
                end2=8) %>% 
  dplyr::mutate(distanceToTSS=end-start_T)

ATAC_peaks_neurons_D868N_HiC_up_dev <- rbind.data.frame(ATAC_peaks_neurons_D868N_HiC_up_dev,ATAC_peaks_neurons_D868N_HiC_up_dev_p)
ATAC_peaks_neurons_D868N_HiC_up_dev  <- ATAC_peaks_neurons_D868N_HiC_up_dev[with(ATAC_peaks_neurons_D868N_HiC_up_dev, ave(IF, peaks, FUN= max)==IF),]


write_delim(ATAC_peaks_neurons_D868N_HiC_up_dev,"SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_UP_HiC_dev",
            delim="\t", col_names = T)

ATAC_peaks_neurons_D868N_HiC_up_dev.bedpe <- ATAC_peaks_neurons_D868N_HiC_up_dev %>% 
  dplyr::select(chr,start,end,chr2,start2,end2,IF) %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_UP_HiC_dev.bedpe",
              delim="\t",col_names = F)

ATAC_peaks_neurons_D868N_HiC_up_dev_symbol <- ATAC_peaks_neurons_D868N_HiC_up_dev %>% # Extracting associated gene list
  dplyr::rename(SYMBOL=8) %>% 
  dplyr::select(SYMBOL) %>% 
  unique %>% 
  write_delim("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_UP_HiC_dev_SYMBOL",
              delim="\t",col_names = T)

Common_coregulated_genes <- inner_join(ATAC_peaks_neurons_D868D_HiC_up_dev_symbol,ATAC_peaks_neurons_D868N_HiC_up_dev_symbol) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/Coregulated_genes_SYMBOL")


ATAC_peaks_neurons_D868D_non_coregulated_SYMBOL <- anti_join(ATAC_peaks_neurons_D868D_HiC_up_dev_symbol,Common_coregulated_genes) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_non_coregulated_SYMBOL")

ATAC_peaks_neurons_D868N_non_coregulated_SYMBOL <- anti_join(ATAC_peaks_neurons_D868N_HiC_up_dev_symbol,Common_coregulated_genes) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_non_coregulated_SYMBOL")


