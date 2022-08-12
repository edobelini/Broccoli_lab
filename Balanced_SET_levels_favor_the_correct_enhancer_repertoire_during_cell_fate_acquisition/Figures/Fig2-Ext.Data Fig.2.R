library(IRanges)
library(dplyr)
library(plyranges)
library(tidyverse)
library(DESeq2)
library(Rsubread)
library(ggpubr)
library(LSD)

# Fig.2 a  Genomic distribution of ATAC cluters peaks


#NPCs D868D


NPC_D868D_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPCD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_D868D_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster1.bed")

NPC_D868D_cluster1_annotate <- ChIPseeker::annotatePeak(NPC_D868D_cluster1,
                                                            tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                            annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster1_annotate.tsv")


table(NPC_D868D_cluster1_annotate$Feature)

NPC_D868D_cluster1_pie <- data.frame(table(NPC_D868D_cluster1_annotate$Feature))

NPC_D868D_cluster1_pie$percentage <- prop.table(NPC_D868D_cluster1_pie$Freq)*100

bp<- ggplot(NPC_D868D_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



NPC_D868D_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPCD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_2"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_D868D_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster2.bed")

NPC_D868D_cluster2_annotate <- ChIPseeker::annotatePeak(NPC_D868D_cluster2,
                                                            tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                            annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster2_annotate.tsv")


NPC_D868D_cluster2_pie <- data.frame(table(NPC_D868D_cluster2_annotate$Feature))

NPC_D868D_cluster2_pie$percentage <- prop.table(NPC_D868D_cluster2_pie$Freq)*100

bp<- ggplot(NPC_D868D_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


NPC_D868D_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPCD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_3"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_D868D_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3.bed")



NPC_D868D_cluster3_annotate <- ChIPseeker::annotatePeak(NPC_D868D_cluster3,
                                                            tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                            annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3_annotate.tsv")

NPC_D868D_cluster3_annotate_SYMBOL <- read_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3_annotate.tsv") %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3_annotate_SYMBOL.tsv")


NPC_D868D_cluster3_pie <- data.frame(table(NPC_D868D_cluster3_annotate$Feature))

NPC_D868D_cluster3_pie$percentage <- prop.table(NPC_D868D_cluster3_pie$Freq)*100

bp<- ggplot(NPC_D868D_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


# Fig.2 e Plot density of cluster 3 peaks inside superenhancers

#load ATAC peaks and Super-enhancers genomic coordinates

cluster3 <- read_bed("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Regions/NPC_D868D_cluster3.bed")

ATAC_NPCD868D_superenhancers <- read_bed("SETBP1_epigenomics/pipeline/Peaks/ATAC_NPCD868D_superenhancers.bed")

Super_enhancers <- read_bed("SETBP1_epigenomics/pipeline/Peaks/NPCD868D_superenhancers.bed")

#Filter out cluster 3 peaks in SE and plot

ATAC_NPCD868D_superenhancers <- pair_overlaps(Super_enhancers,ATAC_NPC_D868D) %>% 
  as.data.frame()

ATAC_NPCD868D_superenhancers$SE <- paste(ATAC_NPCD868D_superenhancers$granges.x.seqnames,ATAC_NPCD868D_superenhancers$granges.x.start,
                                         ATAC_NPCD868D_superenhancers$granges.x.end,sep = ":")

ATAC_NPCD868D_superenhancers_prop <- data.frame(table(ATAC_NPCD868D_superenhancers$SE))

cluster3_Super_enhancers <- pair_overlaps(Super_enhancers,cluster3) %>% 
  as.data.frame() 

cluster3_Super_enhancers$SE <- paste(cluster3_Super_enhancers$granges.x.seqnames,cluster3_Super_enhancers$granges.x.start,
                                     cluster3_Super_enhancers$granges.x.end,sep = ":")

cluster3_Super_enhancers_prop <- data.frame(table(cluster3_Super_enhancers$SE))

ATAC_SE_Comp <- full_join(ATAC_NPCD868D_superenhancers_prop,cluster3_Super_enhancers_prop, by="Var1") %>% 
  replace(., is.na(.), "0") 

ATAC_SE_Comp$Freq.y <- as.numeric(ATAC_SE_Comp$Freq.y)

ATAC_SE_Comp <- ATAC_SE_Comp %>% 
  dplyr::mutate(percentage=(Freq.y/Freq.x)*100) 


ATAC_SE_over_50 <- ATAC_SE_Comp %>% 
  dplyr::filter(percentage>=50) %>% 
  dplyr::rename(SE=1)

write_tsv(ATAC_SE_over_50,"SETBP1_epigenomics/pipeline/Peaks/ATAC_SE_over_50.bed")

ATAC_SE_over_50_peaks <- inner_join(ATAC_SE_over_50,ATAC_NPCD868D_superenhancers, by="SE") %>% 
  dplyr::rename(chr=10,start=11,end=12) %>% 
  dplyr::select(chr,start,end)

write_bed(ATAC_SE_over_50_peaks,"SETBP1_epigenomics/pipeline/Peaks/ATAC_SE_over_50_peaks.bed")

SE_over_50_peaks <- inner_join(ATAC_SE_over_50,ATAC_NPCD868D_superenhancers, by="SE") %>% 
  dplyr::rename(chr=5,start=6,end=7) %>% 
  dplyr::select(chr,start,end) %>% 
  unique() %>% 
  write_bed("SETBP1_epigenomics/pipeline/Peaks/SE_over_50_peaks.bed")

ggplot(ATAC_SE_Comp, aes(x=percentage))+
  geom_density(color="dodgerblue2", fill="cadetblue2")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))
ggsave("SETBP1_epigenomics/pipeline/plots/Cluster3_percentage_in_SE.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 65, units = "mm", dpi = 300, limitsize = TRUE)



# Fig.2 f Quartiles analysis Open Chromatin from PSCs to NPCs

#calculate log2FC and norm counts in Open chromatin of PSCs+ control NPCs using Deseq2

bamsToCount <- dir("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Bowtie2/filtered/PSC_NPC/", full.names = TRUE, pattern = "*.\\.bam$")

consensusToCount <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Regions/PSC_NPC_ATAC_merge.bed",
                             delim="\t",col_names = T) %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) 

consensusToCount

regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount),  
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount,  isPairedEnd = TRUE,
                           countMultiMappingReads = FALSE, maxFragLength = 100,
                           nthreads = 30)

fcResults

myCounts <- fcResults$counts

colnames(myCounts) <- c("NPC_D868D_1","NPC_D868D_2","NPC_D868D_3","NPC_D868N_1","NPC_D868N_2","NPC_D868N_3","PSC_1","PSC_2")



metaData <- data.frame(Group= c("NPC_D868D","NPC_D868D","NPC_D868D","NPC_D868N","NPC_D868N","NPC_D868N","PSC","PSC"),
                       replicates=c("NPC_D868D_1","NPC_D868D_2","NPC_D868D_3","NPC_D868N_1","NPC_D868N_2","NPC_D868N_3","PSC_1","PSC_2"))
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
atacDDS <- DESeq(atacDDS)

NPCMinusH9 <- results(atacDDS, c("Group","NPC_D868D", "PSC"), format = "GRanges")
NPCMinusH9<- NPCMinusH9[order(NPCMinusH9$pvalue)]
  

NPCMinusH9 <- join_overlap_inner(NPCMinusH9,consensusToCount) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  dplyr::rename(GeneID=1)


NPCMinusH9 <- inner_join(NPCMinusH9,atacDDS_normalized_counts) %>% 
  dplyr::mutate(NPC_D868D_counts=(NPC_D868D_1+NPC_D868D_2+NPC_D868D_3)/3,
                NPC_D868N_counts=(NPC_D868N_1+NPC_D868N_2+NPC_D868N_3)/3,
                PSC_counts=(PSC_1+PSC_2)/2)


#Quartiles calculation--------------------------------------------------------------------------------------------------------------------

ATAC_NPC_H9_Q1_only <- NPCMinusH9 %>%
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.25)) %>% 
  dplyr::rename(Q1_NPCvsH9=24) %>%
  dplyr::filter(Q1_NPCvsH9==T) %>% 
  dplyr::rename(chr=2) 

write_bed(ATAC_NPC_H9_Q1_only,"SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q1_only.bed")


ATAC_NPC_H9_Q2_only <- NPCMinusH9 %>%
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.25)) %>% 
  dplyr::rename(Q1_NPCvsH9=24) %>%
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.50))%>%
  dplyr::rename(Q2_NPCvsH9=25) %>% 
  dplyr::filter(Q1_NPCvsH9==F & Q2_NPCvsH9==T) %>% 
  dplyr::rename(chr=2)

write_bed(ATAC_NPC_H9_Q2_only,"SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q2_only.bed")

ATAC_NPC_H9_Q3_only <- NPCMinusH9 %>%
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.25)) %>% 
  dplyr::rename(Q1_NPCvsH9=24) %>%
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.50))%>%
  dplyr::rename(Q2_NPCvsH9=25) %>% 
  dplyr::mutate(log2FoldChange<=quantile(log2FoldChange,0.75))%>%
  dplyr::rename(Q3_NPCvsH9=26)  %>%
  dplyr::filter(Q1_NPCvsH9==F & Q2_NPCvsH9==F & Q3_NPCvsH9==T) %>% 
  dplyr::rename(chr=2)

write_bed(ATAC_NPC_H9_Q3_only,"SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q3_only.bed")

ATAC_NPC_H9_Q4_only <- NPCMinusH9 %>%
  dplyr::mutate(log2FoldChange>quantile(log2FoldChange,0.75))%>%
  dplyr::rename(Q4_NPCvsH9=24) %>%
  dplyr::filter(Q4_NPCvsH9==T) %>% 
  dplyr::rename(chr=2)


write_bed(ATAC_NPC_H9_Q4_only,"SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q4_only.bed")


#Violin plot ATAC in quartiles--------------------------------------------------------------------------

my_comparisons <- list( c("PSCs", "NPCs 
D868D"), c("PSCs", "NPCs 
D868N"), c("NPCs 
D868D", "NPCs 
D868N") )


NPC_PSC_boxplot_Q1 <-ATAC_NPC_H9_Q1_only %>%
  dplyr::select(NPC_D868D_counts,NPC_D868N_counts,PSC_counts) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "PSCs"=3) %>% 
  gather(key=Group, value=Normalized_Counts, "NPCs 
D868D", "NPCs 
D868N", "PSCs")
  
NPC_PSC_boxplot_Q1$Group <- factor(NPC_PSC_boxplot_Q1$Group, levels = c("PSCs","NPCs 
D868D", "NPCs 
D868N"))

ggplot(NPC_PSC_boxplot_Q1) +
  aes(x = Group, y = log2(Normalized_Counts), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#404040","#006d2c","#74c476")) +
  scale_color_manual(values=c("#404040","#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') +
  ylab('log2 (normalized counts)')+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    comparisons = my_comparisons,
    label.y = c(17, 19, 21)  #posizione p value
  ) 

ggsave("SETBP1_epigenomics/pipeline/plots/NPC_PSC_boxplot_Q1.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


NPC_PSC_boxplot_Q2 <-ATAC_NPC_H9_Q2_only  %>%
  dplyr::select(NPC_D868D_counts,NPC_D868N_counts,PSC_counts) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "PSCs"=3) %>% 
  gather(key=Group, value=Normalized_Counts, "NPCs 
D868D", "NPCs 
D868N", "PSCs")

NPC_PSC_boxplot_Q2$Group <- factor(NPC_PSC_boxplot_Q2$Group, levels = c("PSCs","NPCs 
D868D", "NPCs 
D868N"))

ggplot(NPC_PSC_boxplot_Q2) +
  aes(x = Group, y = log2(Normalized_Counts), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#404040","#006d2c","#74c476")) +
  scale_color_manual(values=c("#404040","#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') +
  ylab('log2 (normalized counts)')+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    comparisons = my_comparisons,
    label.y = c(17, 19, 21)  #posizione p value
  ) 

ggsave("SETBP1_epigenomics/pipeline/plots/NPC_PSC_boxplot_Q2.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


NPC_PSC_boxplot_Q3 <-ATAC_NPC_H9_Q3_only %>%
  dplyr::select(NPC_D868D_counts,NPC_D868N_counts,PSC_counts) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "PSCs"=3) %>% 
  gather(key=Group, value=Normalized_Counts, "NPCs 
D868D", "NPCs 
D868N", "PSCs")

NPC_PSC_boxplot_Q3$Group <- factor(NPC_PSC_boxplot_Q3$Group, levels = c("PSCs","NPCs 
D868D", "NPCs 
D868N"))

ggplot(NPC_PSC_boxplot_Q3) +
  aes(x = Group, y = log2(Normalized_Counts), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#404040","#006d2c","#74c476")) +
  scale_color_manual(values=c("#404040","#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') +
  ylab('log2 (normalized counts)')+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    comparisons = my_comparisons,
    label.y = c(17, 19, 21)  #posizione p value
  )  

ggsave("SETBP1_epigenomics/pipeline/plots/NPC_PSC_boxplot_Q3.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


NPC_PSC_boxplot_Q4 <-ATAC_NPC_H9_Q4_only %>%
  dplyr::select(NPC_D868D_counts,NPC_D868N_counts,PSC_counts) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "PSCs"=3) %>% 
  gather(key=Group, value=Normalized_Counts, "NPCs 
D868D", "NPCs 
D868N", "PSCs")

NPC_PSC_boxplot_Q4$Group <- factor(NPC_PSC_boxplot_Q4$Group, levels = c("PSCs","NPCs 
D868D", "NPCs 
D868N"))

ggplot(NPC_PSC_boxplot_Q4) +
  aes(x = Group, y = log2(Normalized_Counts), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#404040","#006d2c","#74c476")) +
  scale_color_manual(values=c("#404040","#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') +
  ylab('log2 (normalized counts)')+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    comparisons = my_comparisons,
    label.y = c(17, 19, 21)  #posizione p value
  ) 

ggsave("SETBP1_epigenomics/pipeline/plots/NPC_PSC_boxplot_Q4.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


#Scatterplot all peaks

h_plot <- NPCMinusH9
h_plot$log_pvalue <- -log(NPCMinusH9$pvalue)

h_plot <- h_plot[!is.na(h_plot$log_pvalue),]

library(plotfunctions)
library(colorRamps)
mycols3 <- matlab.like2(10)
png(filename = "SETBP1_epigenomics/pipeline/plots/All_peaks_PSC_NPC_log2FC.png", width = 25, height = 10, res = 720, units = "cm")
heatscatter(x = h_plot$log2FoldChange, y =h_plot$log_pvalue, cexplot = 0.5, colpal = "matlablike2",main= "", nrcol = 90, grid = 1000, add.contour = F, 
            ylim = c(0, 80), xlim = c(-10, 10), xlab = "Log2 fold change PSCs vs NPCs D868D", ylab = "Log P-value PSCs vs NPCs D868D",frame.plot = FALSE,
            xaxs = "i",
            yaxs = "i")
axis(side=1, at=seq(-10, 10, by=1))
gradientLegend(
  c(-14),
  color = mycols3,
  nCol = 30,
  pos = 0.5,
  side = 4,
  dec = NULL,
  length = 1,
  depth = 0.05,
  inside = FALSE,
  coords = FALSE,
  pos.num = NULL,
  n.seg = NULL,
  border.col = "black",
  tick.col = NULL,
  fit.margin = FALSE
)
legend("topleft", legend=as.vector(SOX2$Motif.Name),  bty = "n", cex=1, pt.cex = 0.7, ncol = 1, y.intersp = 0.5)
dev.off()


#Fig.2 G venn diagram cluster3 peaks in quartiles 

#cluster 3 regions in quartiles----------------------------------------------------------------------------

cluster_3 <- read_bed("SETBP1_epigenomics/pipeline/Peaks/NPC_D868D_cluster3.bed")

cluster_3_top_10k <- read_bed("SETBP1_epigenomics/pipeline/Deseq2/NPC_D868NvsNPC_D868D_down_log2FC_0.bed")

ATAC_NPC_H9_Q1 <- read_bed("SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q1_only.bed")

ATAC_NPC_H9_Q2 <- read_bed("SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q2_only.bed")

ATAC_NPC_H9_Q3 <- read_bed("SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q3_only.bed")

ATAC_NPC_H9_Q4 <- read_bed("SETBP1_epigenomics/pipeline/Peaks/ATAC_NPC_H9_Q4_only.bed")

cluster_3_ATAC_NPC_H9_Q1 <- join_overlap_inner(cluster_3,ATAC_NPC_H9_Q1)

cluster_3_ATAC_NPC_H9_Q2 <- join_overlap_inner(cluster_3,ATAC_NPC_H9_Q2)

cluster_3_ATAC_NPC_H9_Q3 <- join_overlap_inner(cluster_3,ATAC_NPC_H9_Q3)

cluster_3_ATAC_NPC_H9_Q4 <- join_overlap_inner(cluster_3,ATAC_NPC_H9_Q4)

cluster_3_top_10k_ATAC_NPC_H9_Q1 <- join_overlap_inner(cluster_3_top_10k,ATAC_NPC_H9_Q1)

cluster_3_top_10k_ATAC_NPC_H9_Q2 <- join_overlap_inner(cluster_3_top_10k,ATAC_NPC_H9_Q2)

cluster_3_top_10k_ATAC_NPC_H9_Q3 <- join_overlap_inner(cluster_3_top_10k,ATAC_NPC_H9_Q3)

cluster_3_top_10k_ATAC_NPC_H9_Q4 <- join_overlap_inner(cluster_3_top_10k,ATAC_NPC_H9_Q4)

#venn diagram distribution peaks (percentage where calculated at the moment by the operator)

Down_in_NPC <- data.frame(
  group = c("1° quartile", "2° quartile","3° quartile","4° quartile"),
  value = c(4.36,27.48,30.50,37.66)
)
head(df)

bp<- ggplot(Down_in_NPC, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("yellow","orange","red","brown"))+
  scale_color_manual(values=c("yellow","orange","red","brown"))+ theme_void()+
  theme(axis.text.x=element_blank())

ggsave("SETBP1_epigenomics/pipeline/plots/cluster3_in_NPC_D868_quartile_venn.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 80, height = 75, units = "mm", dpi = 300, limitsize = TRUE)  

Down_in_NPC <- data.frame(
  group = c("1° quartile", "2° quartile","3° quartile","4° quartile"),
  value = c(1.62,13.34,26.58,58.46)
)
head(df)

bp<- ggplot(Down_in_NPC, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("yellow","orange","red","brown"))+
  scale_color_manual(values=c("yellow","orange","red","brown"))+ theme_void()+
  theme(axis.text.x=element_blank())

ggsave("SETBP1_epigenomics/pipeline/plots/cluster3_in_NPC_D868_quartile_top10K_venn.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 80, height = 75, units = "mm", dpi = 300, limitsize = TRUE)  



#Extended data Fig.2 a violin plot ATAC in NPCs D868D peaks 


#load Open Chromatin SGS summary

Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T)

#ATAC signal enrichment in NPCs D868D peaks of NPCs D868

Violin_plot_ATAC_NPCs_Ctrl <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(NPC_D868D_peaks==1) %>%
  dplyr::select(NPC_D868D,NPC_D868N) %>% 
  dplyr::rename("NPCs D868D"=1,
                "NPCs D868N"=2) %>% 
  gather(key=Group, value=RPKM, "NPCs D868D","NPCs D868N")

t.test(Violin_plot_ATAC_NPCs_Ctrl)

ggplot(Violin_plot_ATAC_NPCs_Ctrl) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    ref.group = "NPCs D868D",
    label.y = 15  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_NPC_D868_Ctrl_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE) 


#ATAC signal enrichment in NPCs D868D peaks of NPCs I871

Violin_plot_ATAC_NPCs_Ctrl <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(NPC_D868D_peaks==1) %>%
  dplyr::select(NPC_I871I,NPC_I871T) %>% 
  dplyr::rename("NPCs I871I"=1,
                "NPCs I817T"=2) %>% 
  gather(key=Group, value=RPKM, "NPCs I871I","NPCs I817T")

Violin_plot_ATAC_NPCs_Ctrl$Group <- factor(Violin_plot_ATAC_NPCs_Ctrl$Group, levels = c("NPCs I871I","NPCs I817T"))


t.test(Violin_plot_ATAC_NPCs_Ctrl)

ggplot(Violin_plot_ATAC_NPCs_Ctrl) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    ref.group = "NPCs I871I",
    label.y = 15  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_NPC_D868_I871_Ctrl_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)  


#Extended data Fig.2 b c d  Genomic distribution of ATAC cluters peaks

#NPCs SETV5


NPC_SETV5_No_Doxy_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_SETV5_No_Doxy_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster1.bed")

NPC_SETV5_No_Doxy_cluster1_annotate <- ChIPseeker::annotatePeak(NPC_SETV5_No_Doxy_cluster1,
                                                        tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                        annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster1_annotate.tsv")


table(NPC_SETV5_No_Doxy_cluster1_annotate$Feature)

NPC_SETV5_No_Doxy_cluster1_pie <- data.frame(table(NPC_SETV5_No_Doxy_cluster1_annotate$Feature))

NPC_SETV5_No_Doxy_cluster1_pie$percentage <- prop.table(NPC_SETV5_No_Doxy_cluster1_pie$Freq)*100

bp<- ggplot(NPC_SETV5_No_Doxy_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

NPC_SETV5_No_Doxy_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_2"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_SETV5_No_Doxy_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster2.bed")

NPC_SETV5_No_Doxy_cluster2_annotate <- ChIPseeker::annotatePeak(NPC_SETV5_No_Doxy_cluster2,
                                                        tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                        annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster2_annotate.tsv")


NPC_SETV5_No_Doxy_cluster2_pie <- data.frame(table(NPC_SETV5_No_Doxy_cluster2_annotate$Feature))

NPC_SETV5_No_Doxy_cluster2_pie$percentage <- prop.table(NPC_SETV5_No_Doxy_cluster2_pie$Freq)*100

bp<- ggplot(NPC_SETV5_No_Doxy_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


NPC_SETV5_No_Doxy_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NPC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_3"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(NPC_SETV5_No_Doxy_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster3.bed")



NPC_SETV5_No_Doxy_cluster3_annotate <- ChIPseeker::annotatePeak(NPC_SETV5_No_Doxy_cluster3,
                                                        tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                        annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster3_annotate.tsv")

NPC_SETV5_No_Doxy_cluster3_annotate_SYMBOL <- read_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster3_annotate.tsv") %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster3_annotate_SYMBOL.tsv")



NPC_SETV5_No_Doxy_cluster3_pie <- data.frame(table(NPC_SETV5_No_Doxy_cluster3_annotate$Feature))

NPC_SETV5_No_Doxy_cluster3_pie$percentage <- prop.table(NPC_SETV5_No_Doxy_cluster3_pie$Freq)*100

bp<- ggplot(NPC_SETV5_No_Doxy_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_SETV5_No_Doxy_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



#IPSCs SETV5

IPSC_SETV5_No_Doxy_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_IPSC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(IPSC_SETV5_No_Doxy_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster1.bed")

IPSC_SETV5_No_Doxy_cluster1_annotate <- ChIPseeker::annotatePeak(IPSC_SETV5_No_Doxy_cluster1,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                             annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster1_annotate.tsv")


table(IPSC_SETV5_No_Doxy_cluster1_annotate$Feature)

IPSC_SETV5_No_Doxy_cluster1_pie <- data.frame(table(IPSC_SETV5_No_Doxy_cluster1_annotate$Feature))

IPSC_SETV5_No_Doxy_cluster1_pie$percentage <- prop.table(IPSC_SETV5_No_Doxy_cluster1_pie$Freq)*100

bp<- ggplot(IPSC_SETV5_No_Doxy_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


png("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_Common_peaks_cluster1.png", width= 10, height = 10, res=330, units = "cm")
makeVennDiagram(list(IPSC_SETV5_No_Doxy_cluster1, IPSC_SETV5_cluster1), NameOfPeaks=c("", ""),scaled=FALSE, euler.d=FALSE, 
                fill=c("#006d2c", "#74c476"), # circle fill color
                col=c("#006d2c", "#74c476"), #circle border color
                cat.col=c("black", "black"))
dev.off()

IPSC_SETV5_No_Doxy_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_IPSC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_2"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(IPSC_SETV5_No_Doxy_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster2.bed")

IPSC_SETV5_No_Doxy_cluster2_annotate <- ChIPseeker::annotatePeak(IPSC_SETV5_No_Doxy_cluster2,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                             annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster2_annotate.tsv")


IPSC_SETV5_No_Doxy_cluster2_pie <- data.frame(table(IPSC_SETV5_No_Doxy_cluster2_annotate$Feature))

IPSC_SETV5_No_Doxy_cluster2_pie$percentage <- prop.table(IPSC_SETV5_No_Doxy_cluster2_pie$Freq)*100

bp<- ggplot(IPSC_SETV5_No_Doxy_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


IPSC_SETV5_No_Doxy_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_IPSC-SETV5-No_Doxy_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_3"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(IPSC_SETV5_No_Doxy_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster3.bed")



IPSC_SETV5_No_Doxy_cluster3_annotate <- ChIPseeker::annotatePeak(IPSC_SETV5_No_Doxy_cluster3,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                             annoDb="org.Hs.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster3_annotate.tsv")

IPSC_SETV5_No_Doxy_cluster3_annotate_SYMBOL <- read_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster3_annotate.tsv") %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster3_annotate_SYMBOL.tsv")



IPSC_SETV5_No_Doxy_cluster3_pie <- data.frame(table(IPSC_SETV5_No_Doxy_cluster3_annotate$Feature))

IPSC_SETV5_No_Doxy_cluster3_pie$percentage <- prop.table(IPSC_SETV5_No_Doxy_cluster3_pie$Freq)*100

bp<- ggplot(IPSC_SETV5_No_Doxy_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/IPSC_SETV5_No_Doxy_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

#Zebrafish

Zebrafish_GFP_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_Zeb-GFP_ATAC_merge_70_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Zebrafish_GFP_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster1.bed")

Zebrafish_GFP_cluster1_annotate <- ChIPseeker::annotatePeak(Zebrafish_GFP_cluster1,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Drerio.UCSC.danRer11.refGene::TxDb.Drerio.UCSC.danRer11.refGene, 
                                                             annoDb="org.Dr.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster1_annotate.tsv")


table(Zebrafish_GFP_cluster1_annotate$Feature)

Zebrafish_GFP_cluster1_annotate_pie <- data.frame(table(Zebrafish_GFP_cluster1_annotate$Feature))

Zebrafish_GFP_cluster1_annotate_pie$percentage <- prop.table(Zebrafish_GFP_cluster1_annotate_pie$Freq)*100

bp<- ggplot(Zebrafish_GFP_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

Zebrafish_GFP_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_Zeb-GFP_ATAC_merge_70_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Zebrafish_GFP_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster2.bed")

Zebrafish_GFP_cluster1_annotate <- ChIPseeker::annotatePeak(Zebrafish_GFP_cluster1,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Drerio.UCSC.danRer11.refGene::TxDb.Drerio.UCSC.danRer11.refGene, 
                                                             annoDb="org.Dr.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster2_annotate.tsv")


table(Zebrafish_GFP_cluster2_annotate$Feature)

Zebrafish_GFP_cluster2_annotate_pie <- data.frame(table(Zebrafish_GFP_cluster2_annotate$Feature))

Zebrafish_GFP_cluster2_annotate_pie$percentage <- prop.table(Zebrafish_GFP_cluster2_annotate_pie$Freq)*100

bp<- ggplot(Zebrafish_GFP_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



Zebrafish_GFP_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_Zeb-GFP_ATAC_merge_70_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Zebrafish_GFP_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster3.bed")

Zebrafish_GFP_cluster3_annotate <- ChIPseeker::annotatePeak(Zebrafish_GFP_cluster1,
                                                             tssRegion=c(-10000, 2000), TxDb=TxDb.Drerio.UCSC.danRer11.refGene::TxDb.Drerio.UCSC.danRer11.refGene, 
                                                             annoDb="org.Dr.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster1_annotate.tsv")


table(Zebrafish_GFP_cluster3_annotate$Feature)

Zebrafish_GFP_cluster3_annotate_pie <- data.frame(table(Zebrafish_GFP_cluster3_annotate$Feature))

Zebrafish_GFP_cluster3_annotate_pie$percentage <- prop.table(Zebrafish_GFP_cluster3_annotate_pie$Freq)*100

bp<- ggplot(Zebrafish_GFP_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Zebrafish_GFP_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


#Extended Data Fig.2 e Chromatin state annotation NPCs D868D ATAC clusters

NPC_cluster3 <- read_bed("Setbp1_Gdrive/bench/zaghi/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3.bed")



files <- list.files("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Encode_annotation/Hg38/", pattern = "\\.bed", full.names = T)
encode <-  lapply(files, read_delim, delim="\t", skip = 1, col_names = F) 

names(encode) <- str_split_fixed(as.character(list.files("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Encode_annotation/Hg38/", 
                                            pattern = "\\.bed", full.names = F)), pattern = "_", 5)[,1]



over_region <- list()
results <- list()

for (i in 1:length(encode)) {
   join <- join_overlap_inner(NPC_cluster3, encode[[i]] ) %>% 
    as.data.frame()
  over_region[[i]] <- join
  names(over_region)[i] <- names(encode)[i]
  
  join_freq <- data.frame(table(join$X4))
  join_freq$percentage <- prop.table(join_freq$Freq)*100
  join_freq$condition <- paste(names(encode)[i])
  results[[i]] <- join_freq
  names(results)[i] <-  names(encode)[i]
}

prova <- bind_rows(results, .id = "condition") 


prova <- filter(prova, prova$condition %in% c("E007", "E009", "E053", "E081", "E073")) %>% 
  dplyr::rename(Code=4)

prova <- inner_join(prova, metadata)

prova$colour <- c("#FF0000","#C2E105","#C2E105","#FFC34D","#FFC34D","#FFC34D","#FFFF00","#FFFF00","#FFFF00","#FFFF66",
                         "#FF4500","#66CDAA","#FF4500","#8A91D0","#E6B8B7","#7030A0","#808080","gray90", "#FF4500","#FF4500",
                         "#008000", "#008000","#008000","#009600", "#C2E105")

prova$Var1 <- factor(prova$Var1, levels = c("1_TssA","2_PromU","3_PromD1",
                                            "4_PromD2","5_Tx5'","6_Tx",
                                                          "7_Tx3'","8_TxWk","9_TxReg",
                                                          "10_TxEnh5'","11_TxEnh3'","12_TxEnhW",
                                                          "13_EnhA1","14_EnhA2","15_EnhAF",
                                                          "16_EnhW1","17_EnhW2","18_EnhAc",
                                                          "19_DNase","20_ZNF/Rpts","21_Het",
                                                          "22_PromP","23_PromBiv","24_ReprPC",
                                                          "25_Quies"))

bp <- ggplot(prova, aes(x=Descritpion, y=percentage, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")+
  scale_fill_manual( values=c("#FF0000","#66CDAA","#FF4500",
                       "#FF4500","#008000","#008000",
                       "#008000","#009600","#C2E105",
                       "#C2E105","#C2E105","#FFC34D",
                       "#FFC34D","#FFC34D","#FFFF00",
                       "#FFFF00","#FFFF00","#FFFF66",
                       "#FF4500","#FF4500","#8A91D0",
                       "#E6B8B7","#7030A0","#808080",
                       "gray90"))+
  scale_color_manual( values=c("#FF0000","#66CDAA","#FF4500",
                        "#FF4500","#008000","#008000",
                        "#008000","#009600","#C2E105",
                        "#C2E105","#C2E105","#FFC34D",
                        "#FFC34D","#FFC34D","#FFFF00",
                        "#FFFF00","#FFFF00","#FFFF66",
                        "#FF4500","#FF4500","#8A91D0",
                        "#E6B8B7","#7030A0","#808080",
                        "gray90")) +
  xlab("")+
  ylab("percentage")+
  coord_flip()+
  guides(fill=guide_legend(title="chromHMM state")) +
  theme_classic()+ theme(axis.text.x = element_text(size = 34,family = "Arial"),
          axis.text.y = element_text(size = 34,family = "Arial"),
          axis.title.x = element_text(size = 34,family = "Arial"),
          axis.line = element_line(size = 2),
          legend.text=element_text(size = 34,family = "Arial"),
          legend.title = element_text(size = 34,family = "Arial")) + theme(legend.position = "none")

bp 

ggsave("SETBP1_epigenomics/pipeline/plots/Cluster_3_encode.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 700, height = 305, units = "mm", dpi = 300, limitsize = TRUE) 

wgEncodeRegTfbsClustered <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Encode_annotation/Hg38/wgEncodeRegTfbsClusteredV3_hg38.bed") %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                CHIP=4)

CHIP <- as.data.frame(table(wgEncodeRegTfbsClustered$CHIP))

CHIP <- as.character(CHIP$Var1)

NPC_cluster3 <- read_bed("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster3.bed")

pt_res <- list()

ENCODE <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Encode_annotation/Hg38/wgEncodeRegTfbsClusteredV3_hg38.bed") %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                CHIP=4)

CHIP <- as.data.frame(table(ENCODE$CHIP))

CHIP <- as.character(CHIP$Var1)

genome <- filterChromosomes(getGenome("hg38")) 
genome=genome

for (i in 1:length(CHIP)) {
  ENCODE_F <- ENCODE %>% 
    dplyr::filter(CHIP==c((paste("",CHIP[i],sep = ""))))%>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
pt <- overlapPermTest(A=NPC_cluster3, B=ENCODE_F, force.parallel=T, min.parallel=10,ntimes=10000, genome=genome)
  pt_res[[i]] <- data.frame(CHIP_NAME=paste("",CHIP[i],sep = ""),
                            z_score=pt$numOverlaps$zscore,
                            p_val=pt$numOverlaps$pval)
  names(pt_res)[i] <- paste("",CHIP[i],sep = "")
  print(paste("Making permutation of sample", "",CHIP[i],sep = ""))
}


pt_res_df <- bind_rows(pt_res, .id = "CHIP_NAME") 
plot(pt_res_df$z_score, pt_res_df$p_val)

NPC_cluster2 <- read_bed("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster2.bed")

over_region <- list()
results <- list()

for (i in 1:length(encode)) {
  join <- join_overlap_inner(NPC_cluster2, encode[[i]]) %>% 
    as.data.frame()
  over_region[[i]] <- join
  names(over_region)[i] <- names(encode)[i]
  
  join_freq <- data.frame(table(join$X4))
  join_freq$percentage <- prop.table(join_freq$Freq)*100
  join_freq$condition <- paste(names(encode)[i])
  results[[i]] <- join_freq
  names(results)[i] <-  names(encode)[i]
}

prova <- bind_rows(results, .id = "condition") 


prova <- filter(prova, prova$condition %in% c("E007", "E009", "E053", "E081", "E073")) %>% 
  dplyr::rename(Code=4)

prova <- inner_join(prova, metadata)

prova$Var1 <- factor(prova$Var1, levels = c("1_TssA","2_PromU","3_PromD1",
                                            "4_PromD2","5_Tx5'","6_Tx",
                                            "7_Tx3'","8_TxWk","9_TxReg",
                                            "10_TxEnh5'","11_TxEnh3'","12_TxEnhW",
                                            "13_EnhA1","14_EnhA2","15_EnhAF",
                                            "16_EnhW1","17_EnhW2","18_EnhAc",
                                            "19_DNase","20_ZNF/Rpts","21_Het",
                                            "22_PromP","23_PromBiv","24_ReprPC",
                                            "25_Quies"))


bp <- ggplot(prova, aes(x=Descritpion, y=percentage, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")+
  scale_fill_manual(values=c("#FF0000","#66CDAA","#FF4500",
                             "#FF4500","#008000","#008000",
                             "#008000","#009600","#C2E105",
                             "#C2E105","#C2E105","#FFC34D",
                             "#FFC34D","#FFC34D","#FFFF00",
                             "#FFFF00","#FFFF00","#FFFF66",
                             "#FF4500","#FF4500","#8A91D0",
                             "#E6B8B7","#7030A0","#808080",
                             "gray90"))+
  scale_color_manual(values=c("#FF0000","#66CDAA","#FF4500",
                              "#FF4500","#008000","#008000",
                              "#008000","#009600","#C2E105",
                              "#C2E105","#C2E105","#FFC34D",
                              "#FFC34D","#FFC34D","#FFFF00",
                              "#FFFF00","#FFFF00","#FFFF66",
                              "#FF4500","#FF4500","#8A91D0",
                              "#E6B8B7","#7030A0","#808080",
                              "gray90")) +
  xlab("")+
  ylab("percentage")+
  coord_flip()+
  guides(fill=guide_legend(title="chromHMM state")) +
  theme_classic()+ theme(axis.text.x = element_text(size = 34,family = "Arial"),
                         axis.text.y = element_text(size = 34,family = "Arial"),
                         axis.title.x = element_text(size = 34,family = "Arial"),
                         axis.line = element_line(size = 2),
                         legend.text=element_text(size = 34,family = "Arial"),
                         legend.title = element_text(size = 34,family = "Arial")) + theme(legend.position = "none")
bp 

ggsave("SETBP1_epigenomics/pipeline/plots/Cluster_2_encode.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 700, height = 305, units = "mm", dpi = 300, limitsize = TRUE) 

NPC_cluster1 <- read_bed("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/NPC_D868D_cluster1.bed")

over_region <- list()
results <- list()


for (i in 1:length(encode)) {
  join <- join_overlap_inner(NPC_cluster1, encode[[i]]) %>% 
    as.data.frame()
  over_region[[i]] <- join
  names(over_region)[i] <- names(encode)[i]
  
  join_freq <- data.frame(table(join$X4))
  join_freq$percentage <- prop.table(join_freq$Freq)*100
  join_freq$condition <- paste(names(encode)[i])
  results[[i]] <- join_freq
  names(results)[i] <-  names(encode)[i]
}

prova <- bind_rows(results, .id = "condition") 


prova <- filter(prova, prova$condition %in% c("E007", "E009", "E053", "E081", "E073")) %>% 
  dplyr::rename(Code=4)

prova <- inner_join(prova, metadata)

prova$Var1 <- factor(prova$Var1, levels = c("1_TssA","2_PromU","3_PromD1",
                                            "4_PromD2","5_Tx5'","6_Tx",
                                            "7_Tx3'","8_TxWk","9_TxReg",
                                            "10_TxEnh5'","11_TxEnh3'","12_TxEnhW",
                                            "13_EnhA1","14_EnhA2","15_EnhAF",
                                            "16_EnhW1","17_EnhW2","18_EnhAc",
                                            "19_DNase","20_ZNF/Rpts","21_Het",
                                            "22_PromP","23_PromBiv","24_ReprPC",
                                            "25_Quies"))

bp <- ggplot(prova, aes(x=Descritpion, y=percentage, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")+
  scale_fill_manual(values=c("#FF0000","#66CDAA","#FF4500",
                             "#FF4500","#008000","#008000",
                             "#008000","#009600","#C2E105",
                             "#C2E105","#C2E105","#FFC34D",
                             "#FFC34D","#FFC34D","#FFFF00",
                             "#FFFF00","#FFFF00","#FFFF66",
                             "#FF4500","#FF4500","#8A91D0",
                             "#E6B8B7","#7030A0","#808080",
                             "gray90"))+
  scale_color_manual(values=c("#FF0000","#66CDAA","#FF4500",
                              "#FF4500","#008000","#008000",
                              "#008000","#009600","#C2E105",
                              "#C2E105","#C2E105","#FFC34D",
                              "#FFC34D","#FFC34D","#FFFF00",
                              "#FFFF00","#FFFF00","#FFFF66",
                              "#FF4500","#FF4500","#8A91D0",
                              "#E6B8B7","#7030A0","#808080",
                              "gray90"))+
  xlab("")+
  ylab("percentage")+
  coord_flip()+
  guides(fill=guide_legend(title="chromHMM state")) +
  theme_classic()+ theme(axis.text.x = element_text(size = 34,family = "Arial"),
                         axis.text.y = element_text(size = 34,family = "Arial"),
                         axis.title.x = element_text(size = 34,family = "Arial"),
                         axis.line = element_line(size = 2),
                         legend.text=element_text(size = 34,family = "Arial"),
                         legend.title = element_text(size = 34,family = "Arial")) + theme(legend.position = "none")

bp 

ggsave("SETBP1_epigenomics/pipeline/plots/Cluster_1_encode.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 700, height = 305, units = "mm", dpi = 300, limitsize = TRUE) 

