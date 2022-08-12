library(IRanges)
library(dplyr)
library(plyranges)
library(tidyverse)





# Fig.2 Plot density of cluster 3 peaks inside superenhancers

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

       