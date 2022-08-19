library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(IRanges)
library(dplyr)
library(plyranges)
library(remotes)
library(DESeq2)
library(Rsubread)
library(Tidyverse)
library(dplyr)



#Bam files list where to perform Counts calculation

bamsToCount <- dir("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/Neuro_dev/", full.names = TRUE, pattern = "*.\\.bam$")

# Calculate differential accessibility using DESEQ2

consensusToCount <- read_bed("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Regions/SETBP1_ATAC_merge.bed") #list of regions where to perform calculation


regionsToCount <- data.frame(GeneID = paste(seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), chr = seqnames(consensusToCount), 
                             start = start(consensusToCount), end = end(consensusToCount), Strand = strand(consensusToCount)) 


fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE,
                           countMultiMappingReads = FALSE, maxFragLength = 100,
                           nthreads = 30)


myCounts <- fcResults$counts

colnames(myCounts) <- c("NPC_D868D_1","NPC_D868D_2","NPC_D868D_3",
                        "NPC_D868N_1","NPC_D868N_2","NPC_D868N_3",
                        "Neu_D868D_1","Neu_D868D_2","Neu_D868D_3",
                        "Neu_D868N_1","Neu_D868N_2","Neu_D868N_3")



metaData<- data.frame(Group= c("NPC_D868D","NPC_D868D","NPC_D868D",
                               "NPC_D868N","NPC_D868N","NPC_D868N",
                               "Neu_D868D","Neu_D868D","Neu_D868D",
                               "Neu_D868N","Neu_D868N","Neu_D868N"),
                                replicates=c("NPC_D868D_1","NPC_D868D_2","NPC_D868D_3",
                                             "NPC_D868N_1","NPC_D868N_2","NPC_D868N_3",
                                             "Neu_D868D_1","Neu_D868D_2","Neu_D868D_3",
                                             "Neu_D868N_1","Neu_D868N_2","Neu_D868N_3"))

atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, design = ~Group, rowRanges = consensusToCount)

atacDDS <- DESeq(atacDDS)
atacDDS_normalized_counts <- counts(atacDDS,normalized=T) %>% 
  as.data.frame()

atacDDS_normalized_counts <- rownames_to_column(atacDDS_normalized_counts) %>% 
  dplyr::rename(GeneID=1)


atacDDS_normalized_counts <- inner_join(atacDDS_normalized_counts,regionsToCount) %>% 
  write_delim("SETBP1_epigenomics/pipeline/Deseq2/atacDDS_normalized_counts",
              delim="\t",col_names = T)


#calculate differential peaks 

#Neu_D868NvsNeu_D868D comparison

Neu_D868DvsNPC_D868D <- results(atacDDS, c("Group", "Neu_D868D", "NPC_D868D"), format = "GRanges")
Neu_D868DvsNPC_D868D<- Neu_D868DvsNPC_D868D[order(Neu_D868DvsNPC_D868D$pvalue)] %>% 
  as.data.frame()

Neu_D868DvsNPC_D868D_up <- Neu_D868DvsNPC_D868D %>% 
  dplyr::filter(log2FoldChange> 0.5 & padj <= 0.05)%>% 
  dplyr::rename(chr=1)

write_tsv(Neu_D868DvsNPC_D868D_up,"SETBP1_epigenomics/pipeline/Deseq2/Neu_D868DvsNPC_D868D_UP.tsv")

write_bed(Neu_D868DvsNPC_D868D_up,"SETBP1_epigenomics/pipeline/Deseq2/Neu_D868DvsNPC_D868D_UP.bed")


Neu_D868DvsNPC_D868D_down <- Neu_D868DvsNPC_D868D %>% 
  dplyr::filter(log2FoldChange< -0.5 & padj <= 0.05)%>% 
  dplyr::rename(chr=1)
write_tsv(Neu_D868DvsNPC_D868D_down,"SETBP1_epigenomics/pipeline/Deseq2/Neu_D868DvsNPC_D868D_down.tsv")

write_bed(Neu_D868DvsNPC_D868D_down,"SETBP1_epigenomics/pipeline/Deseq2/Neu_D868DvsNPC_D868D_down.bed")