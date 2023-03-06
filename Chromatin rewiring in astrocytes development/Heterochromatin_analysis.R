library(tidyverse)
library(esquisse)
library(clipr)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(IRanges)
library(dplyr)
library(plyranges)
library(remotes)
library(DESeq2)
library(Rsubread)
library(ComplexHeatmap)
library(circlize)
library(gtools)
library(EnhancedVolcano)
library(ggplot2)
library(janitor)
library(future)
library(ggpubr)
library(ChIPseeker)
library(rGREAT)
library(ChIPpeakAnno)
library(data.table)
library(readxl)
library(openxlsx)
library(gprofiler2)
library(ggnewscale)
setwd("Cut_Tag_Astro/R_analysis/")


#Isolate cluster 3 of H3K27me3 regions which show an increse in methylation in mature astrocytes

H3K27me3_2M <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/Cut_Tag_Astro/Peak_centered/Astro_all_in_Astro_H3K27me3_IGV_40_median_3_Compute_Matrix_heatmap.bed") %>% 
  dplyr::rename(chr=1) %>% 
  dplyr::filter(deepTools_group==c("cluster_3")) %>% 
  write_bed("Cut_Tag_Astro/R_analysis/cluster3_H3K27me3.bed")

cluster3_H3K27me3_annotate <- ChIPseeker::annotatePeak("Cut_Tag_Astro/R_analysis/cluster3_H3K27me3",
                                                      tssRegion=c(-3000, 3000), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                      annoDb="org.Mm.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>%
  write_tsv("Cut_Tag_Astro/R_analysis/cluster3_H3K27me3_annotate.tsv")

cluster3_H3K27me3_annotate_pie <- data.frame(table(cluster3_H3K27me3_annotate$Feature))

cluster3_H3K27me3_annotate_pie$percentage <- prop.table(cluster3_H3K27me3_annotate_pie$Freq)*100

bp<- ggplot(cluster3_H3K27me3_annotate, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868N_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



library(gprofiler2)

gostres <- gost(query = cluster3_H3K27me3_annotate$SYMBOL,
                organism = "mmusculus",
                evcodes = TRUE,
                significant = TRUE,
                correction_method = "fdr",
                user_threshold = 0.05 , sources = c("GO:BP"))

GO <- as.data.frame(gostres$result)
GO$`Percentage of enrichment` <- GO$intersection_size / GO$term_size *100
GO$`-log10 Pvalue` <- -log10(GO$p_value) 

write.xlsx(GO,"Cut_Tag_Astro/R_analysis/GO_cluster3_H3K27me3.xlsx")

GO$Condition <- "H3K27me3 Cluster3"

GO_top10 <- GO[c(1:10),]

GO_top10 %>%
  ggplot() +
  scale_size(range = c(1, 4))+
  geom_point(
    aes(x=Condition, y=reorder(term_name, `-log10 Pvalue`), size=`Percentage of enrichment`, color=`-log10 Pvalue`))+
  #scale_color_gradientn(colors = rev(brewer.ylorrd(100)))+
  scale_color_gradient2(low ="grey70", high = "#006d2c", mid = "#CCFF99", midpoint = -1.8) +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=10, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 0),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=10),
        axis.title.y = element_text(face = "bold", color = "black", size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        legend.position="left",
        panel.background = element_rect(fill = "white",colour = "white", size = 1, linetype = "solid"))+
  ggtitle("GO Biological processes")

GO_neu <- GO[GO$term_name %like% c("neu"),]

GO_neu_top10 <- GO_neu[c(1:10),]

GO_neu_top10 %>%
  ggplot() +
  scale_size(range = c(1, 4))+
  geom_point(
    aes(x=Condition, y=reorder(term_name, `-log10 Pvalue`), size=`Percentage of enrichment`, color=`-log10 Pvalue`))+
  #scale_color_gradientn(colors = rev(brewer.ylorrd(100)))+
  scale_color_gradient2(low ="grey70", high = "#006d2c", mid = "#CCFF99", midpoint = -1.8) +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=10, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 0),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=10),
        axis.title.y = element_text(face = "bold", color = "black", size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        legend.position="left",
        panel.background = element_rect(fill = "white",colour = "white", size = 1, linetype = "solid"))+
  ggtitle("GO Biological processes")

#Isolate cluster 3 of H3K9me3 regions which show an increse in methylation in mature astrocytes

H3K9me3_2M <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/Cut_Tag_Astro/Peak_centered/Astro_all_in_Astro_H3K9me3_IGV_40_median_3_Compute_Matrix_heatmap.bed") %>% 
  dplyr::rename(chr=1) %>% 
  dplyr::filter(deepTools_group==c("cluster_3")) %>% 
  write_tsv("Cut_Tag_Astro/R_analysis/cluster3_H3K9me3")

H3K9me3_2M <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/Cut_Tag_Astro/Peak_centered/Astro_all_in_Astro_H3K9me3_IGV_40_median_3_Compute_Matrix_heatmap.bed") %>% 
  dplyr::rename(chr=1) %>% 
  dplyr::filter(deepTools_group==c("cluster_3")) %>% 
  write_bed("Cut_Tag_Astro/R_analysis/cluster3_H3K9me3.bed")

cluster3_H3K9me3_annotate <- ChIPseeker::annotatePeak("Cut_Tag_Astro/R_analysis/cluster3_H3K9me3",
                                                      tssRegion=c(-3000, 3000), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                      annoDb="org.Mm.eg.db") %>% 
  as_tibble() %>% 
  dplyr::rename(chr=seqnames,ensembl_gene_id=ENSEMBL) %>% 
  mutate(Feature=annotation,
         Feature=gsub(" \\(.*","",Feature),
         Feature=gsub("Distal Intergenic","Intergenic",Feature),
         Feature=gsub("3' UTR","Exon",Feature),
         Feature=gsub("5' UTR","Exon",Feature),
         Feature=gsub("Downstream","Intergenic",Feature)) %>% 
  replace(., is.na(.), "0") %>%
  write_tsv("Cut_Tag_Astro/R_analysis/cluster3_H3K9me3_annotate.tsv")

cluster3_H3K9me3_annotate_pie <- data.frame(table(cluster3_H3K9me3_annotate$Feature))

cluster3_H3K9me3_annotate_pie$percentage <- prop.table(cluster3_H3K9me3_annotate_pie$Freq)*100

bp<- ggplot(cluster3_H3K9me3_annotate_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Cut_Tag_Astro/plots/cluster3_H3K9me3_annotate_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



library(gprofiler2)

gostres <- gost(query = cluster3_H3K9me3_annotate$SYMBOL,
                organism = "mmusculus",
                evcodes = TRUE,
                significant = TRUE,
                correction_method = "fdr",
                user_threshold = 0.05 , sources = c("GO:BP"))

GO <- as.data.frame(gostres$result)
GO$`Percentage of enrichment` <- GO$intersection_size / GO$term_size *100
GO$`-log10 Pvalue` <- -log10(GO$p_value) 

write.xlsx(GO,"Cut_Tag_Astro/R_analysis/GO_cluster3_H3K9me3.xlsx")

GO$Condition <- "H3K9me3 Cluster3"

GO_neu <- GO[GO$term_name %like% c("neu"),]

GO_neu_top20 <- GO_neu[c(1:10),]

GO_neu_top20 %>%
  ggplot() +
  scale_size(range = c(1, 4))+
  geom_point(
    aes(x=Condition, y=reorder(term_name, `-log10 Pvalue`), size=`Percentage of enrichment`, color=`-log10 Pvalue`))+
  #scale_color_gradientn(colors = rev(brewer.ylorrd(100)))+
  scale_color_gradient2(low ="grey70", high = "#08519c", mid = "#99CCFF", midpoint = -1.8) +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=10, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 0),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=10),
        axis.title.y = element_text(face = "bold", color = "black", size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        legend.position="left",
        panel.background = element_rect(fill = "white",colour = "white", size = 1, linetype = "solid"))+
  ggtitle("GO Biological processes")



#load RNA-seq data to perform interpolation between RNA-seq and cluster3

degs_astro <- read_xlsx("Share_HSR/Ric.Broccoli/zaghi.mattia/Cut_Tag_Astro/RNA_seq/dge_normcounts_astro.xlsx") %>% 
  dplyr::rename(P_1=1,P_2=2,P_3=3,
                M_1=4,M_2=5,M_3=6,
                SYMBOL=7) %>% 
  dplyr::mutate(P=(P_1+P_2+P_3)/3,
                M=(M_1+M_2+M_3)/3) %>% 
  dplyr::filter(padj<=0.05)
  

degs_astro_up <- degs_astro %>% 
  dplyr::filter(log2FoldChange>0)

degs_astro_down <- degs_astro %>% 
  dplyr::filter(log2FoldChange<0)



