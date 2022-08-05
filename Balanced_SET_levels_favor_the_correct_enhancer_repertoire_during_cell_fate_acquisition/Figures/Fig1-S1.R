library(dplyr)
library(tidyverse)
library(ggpubr)
library(clipr)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(IRanges)
library(dplyr)
library(plyranges)
library(ChIPseeker)
library(org.Hs.eg.db)

#Load multiBigWigSummary containig all annotation relative to H3K27ac

multiBigWigSummary <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_Chip_H3K27ac_OpenChromatin_table_annotate", col_names = F) %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                NPC_D868D_ATAC=4,
                NPC_D868N_ATAC=5,
                NPC_D868D_H3K27ac=6,
                NPC_D868N_H3K27ac=7,
                NPC_D868D_SET=8,
                NPC_D868N_SET=9,
                NPC_D868D=10,
                NPC_D868N=11) %>% 
  dplyr::mutate(NPC_D868D_H3K27ac_peaks=NPC_D868D_H3K27ac_peaks/NPC_D868D_H3K27ac_peaks,
                NPC_D868N_H3K27ac_peaks=NPC_D868N_H3K27ac_peaks/NPC_D868N_H3K27ac_peaks) %>% 
  mutate_each(funs(replace(., is.na(.), 0))) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_Chip_H3K27ac_OpenChromatin_table_annotate")


#Load multiBigWigSummary containig all annotation of SGS ATAC peaks patients
Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                NPC_D868D=4,
                NPC_D868N=5,
                NPC_I871I=6,
                NPC_I871T=7,
                Neu_D868D=8,
                Neu_D868N=9,
                dir_Neu_D868D=10,
                dir_Neu_D868N=11,
                NPC_D868D_peaks=12,
                NPC_D868N_peaks=13,
                NPC_I871I_peaks=14,
                NPC_I871T_peaks=15,
                Neu_D868D_peaks=16,
                Neu_D868N_peaks=17,
                dir_Neu_D868D_peaks=18,
                dir_Neu_D868N_peaks=19,
                Neu_D868D_up_dev=20,
                Neu_D868D_down_dev=21,
                Neu_D868N_up_dev=22,
                Neu_D868N_down_dev=23) %>% 
  dplyr::mutate( NPC_D868D_peaks=NPC_D868D_peaks/NPC_D868D_peaks,
                 NPC_D868N_peaks=NPC_D868N_peaks/NPC_D868N_peaks,
                 NPC_I871I_peaks=NPC_I871I_peaks/NPC_I871I_peaks,
                 NPC_I871T_peaks=NPC_I871T_peaks/NPC_I871T_peaks,
                 Neu_D868D_peaks=Neu_D868D_peaks/Neu_D868D_peaks,
                 Neu_D868N_peaks=Neu_D868N_peaks/Neu_D868N_peaks,
                 dir_Neu_D868D_peaks=dir_Neu_D868D_peaks/dir_Neu_D868D_peaks,
                 dir_Neu_D868N_peaks=dir_Neu_D868N_peaks/dir_Neu_D868N_peaks) %>% 
  mutate_each(funs(replace(., is.na(.), 0))) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate")


#Fig.1 d Violin plot H3K27ac in all NPC D868D H3K27ac peaks

Box_plot_NPC_D868D_K27ac_ctrl <- multiBigWigSummary %>% 
  dplyr::filter(NPC_D868D_H3K27ac_peaks==1) %>% 
  dplyr::select(NPC_D868D_H3K27ac,NPC_D868N_H3K27ac) %>% 
  dplyr::rename("NPC D868D"=1,
                "NPC D868N"=2) %>% 
  gather(key=Group, value=RPKM, "NPC D868D","NPC D868N")

ggplot(Box_plot_NPC_D868D_K27ac_ctrl) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') +
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"),
        axis.line = element_line(size = 1)) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    ref.group = "NPC D868D",
    label.y = 15  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Box_plot_NPC_H3K27ac_ctrl.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)

#Fig.1 e Violin plot ATAC in all NPC D868D H3K27ac peaks

Violin_NPC_Ctrl_ATAC_in_H3K27ac <- multiBigWigSummary %>%
  dplyr::filter(NPC_D868D_H3K27ac_peaks==1) %>% 
  dplyr::select(NPC_D868D_ATAC, NPC_D868N_ATAC) %>% 
  dplyr::rename("NPC D868D"=1,
                "NPC D868N"=2) %>% 
  gather(key=Group, value=RPKM, "NPC D868D","NPC D868N")

ggplot(Violin_NPC_Ctrl_ATAC_in_H3K27ac) +
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
               colour = "black", width = 0.05) +
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    ref.group = "NPC D868D",
    label.y = 16  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_NPC_Ctrl_ATAC_in_H3K27ac.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)  


#Fig.1 g Violin plot ATAC in all NPC D868 ATAC peaks

Violin_plot_ATAC_NPCs_all <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(NPC_D868D_peaks==1 | NPC_D868N_peaks==1) %>%
  dplyr::select(NPC_D868D,NPC_D868N) %>% 
  dplyr::rename("NPCs D868D"=1,
                "NPCs D868N"=2) %>% 
  gather(key=Group, value=RPKM, "NPCs D868D","NPCs D868N")

ggplot(Violin_plot_ATAC_NPCs_all) +
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

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_NPC_D868_all_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE) 


#Fig.1 i Violin plot ATAC in all NPC I871 peaks 

Violin_plot_ATAC_NPCs_all <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(NPC_I871I_peaks==1 | NPC_I871T_peaks==1) %>%
  dplyr::select(NPC_I871I,NPC_I871T) %>% 
  dplyr::rename("NPC I871I"=1,
                "NPC I871T"=2) %>% 
  gather(key=Group, value=RPKM, "NPC I871I","NPC I871T") 

ggplot(Violin_plot_ATAC_NPCs_all) +
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
    ref.group = "NPC I871I",
    label.y = 16  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_NPC_I871_all_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)


#Fig.S1 h Correlation plot between ATAC signal, RNA and H3K27ac in NPC D868D & NPC D868N

#ATAC-RNA Correlation

RNA_seq <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/RNA_seq/Results_pipeline_hg38/normalized_counts.tsv") %>% 
  dplyr::select(Gene,NPCs_D868D_1,NPCs_D868D_2,NPCs_D868D_6,NPCs_D868N_5, NPCs_D868N_7, NPCs_D868N_8) %>% 
  dplyr::mutate(NPC_D868D_RNA=(NPCs_D868D_1+NPCs_D868D_2+NPCs_D868D_6)/3,
                NPC_D868N_RNA=(NPCs_D868N_5+NPCs_D868N_7+NPCs_D868N_8)/3) %>% 
  dplyr::rename(SYMBOL=1)


NPC_D868D_annotate <- ChIPseeker::annotatePeak("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate",
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
  dplyr::filter(NPC_D868D_peaks==1) %>% 
  dplyr::select(NPC_D868D,SYMBOL)

NPC_D868D_annotate_RNA_seq <- inner_join(NPC_D868D_annotate,
                                         RNA_seq) %>% 
  dplyr::select(NPC_D868D,NPC_D868D_RNA,SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  summarise(NPC_D868D=log2(mean(NPC_D868D)),NPC_D868D_RNA=log2(mean(NPC_D868D_RNA)))
  

ggplot(NPC_D868D_annotate_RNA_seq, aes(x = NPC_D868D, y = NPC_D868D_RNA)) + 
  geom_point() + 
  geom_abline(colour = "blue")+
  stat_cor(method = "pearson", label.x = 7, label.y = 22, size=6)+
  xlab('ATAC log2(RPKM) NPCs D868D') +
  ylab('RNA log2(nor.counts) NPCs D868D')+
  theme_classic()+
  xlim(0,15)+
  theme(axis.text.y = element_text(size = 30,family = "Arial"),
        axis.title.y = element_text(size = 28,family = "Arial"),
        axis.title.x = element_text(size = 28,family = "Arial"),
        axis.text.x = element_text(size = 30,family = "Arial"),
        axis.line = element_line(size = 2))

ggsave("SETBP1_epigenomics/pipeline/plots/Corr_ATAC_RNA_NPC_D868D.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 300, height = 250, units = "mm", dpi = 300, limitsize = TRUE)  


NPC_D868N_annotate <- ChIPseeker::annotatePeak("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate",
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
  dplyr::filter(NPC_D868N_peaks==1) %>% 
  dplyr::select(NPC_D868N,SYMBOL)

NPC_D868N_annotate_RNA_seq <- inner_join(NPC_D868N_annotate,
                                         RNA_seq) %>% 
  dplyr::select(NPC_D868N,NPC_D868N_RNA,SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  summarise(NPC_D868N=log2(mean(NPC_D868N)),NPC_D868N_RNA=log2(mean(NPC_D868N_RNA)))


ggplot(NPC_D868N_annotate_RNA_seq, aes(x = NPC_D868N, y = NPC_D868N_RNA)) + 
  geom_point() + 
  geom_abline(colour = "blue")+
  stat_cor(method = "pearson", label.x = 7, label.y = 22, size=6)+
  xlab('ATAC log2(RPKM) NPCs D868N') +
  ylab('RNA log2(nor.counts) NPCs D868N')+
  theme_classic()+
  xlim(0,15)+
  theme(axis.text.y = element_text(size = 30,family = "Arial"),
        axis.title.y = element_text(size = 28,family = "Arial"),
        axis.title.x = element_text(size = 28,family = "Arial"),
        axis.text.x = element_text(size = 30,family = "Arial"),
        axis.line = element_line(size = 2))


ggsave("SETBP1_epigenomics/pipeline/plots/Corr_ATAC_RNA_NPC_D868N.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 300, height = 250, units = "mm", dpi = 300, limitsize = TRUE)


#ATAC-H3K27ac Correlation
ATAC_H3K27ac_NPC_D868 <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_H3K27ac_NPC_D868_merge_annotate.bed", col_names = T) %>% 
  dplyr::rename(chr=1,
                start=2,
                end=3,
                NPC_D868D=4,
                NPC_D868N=5,
                NPC_D868D_H3K27ac=6,
                NPC_D868N_H3K27ac=7,
                NPC_D868D_ATAC_peaks=8,
                NPC_D868N_ATAC_peaks=9,
                NPC_D868D_H3K27ac_peaks=10,
                NPC_D868N_H3K27ac_peaks=11) %>% 
  dplyr::mutate(NPC_D868D_ATAC_peaks=NPC_D868D_ATAC_peaks/NPC_D868D_ATAC_peaks,
                NPC_D868N_ATAC_peaks=NPC_D868N_ATAC_peaks/NPC_D868N_ATAC_peaks,
                NPC_D868D_H3K27ac_peaks=NPC_D868D_H3K27ac_peaks/NPC_D868D_H3K27ac_peaks,
                NPC_D868N_H3K27ac_peaks=NPC_D868N_H3K27ac_peaks/NPC_D868N_H3K27ac_peaks) %>% 
  mutate_each(funs(replace(., is.na(.), 0)))

ATAC_H3K27ac_NPC_D868D <- ATAC_H3K27ac_NPC_D868  %>% 
  dplyr::filter(NPC_D868D_ATAC_peaks==1 & NPC_D868D_H3K27ac_peaks==1) %>% 
  dplyr::select(NPC_D868D,NPC_D868D_H3K27ac)

ggplot(ATAC_H3K27ac_NPC_D868D, aes(x = log2(NPC_D868D), y = log2(NPC_D868D_H3K27ac))) + 
  geom_point() + 
  geom_abline(colour = "blue")+
  stat_cor(method = "pearson", label.x = 7, label.y = 13, size=6)+
  xlab('ATAC log2(RPKM) NPCs D868D') +
  ylab('H3K27ac log2(RPKM) NPCs D868D')+
  theme_classic()+
  ylim(0,15)+
  xlim(0,13)+
  theme(axis.text.y = element_text(size = 30,family = "Arial"),
        axis.title.y = element_text(size = 28,family = "Arial"),
        axis.title.x = element_text(size = 28,family = "Arial"),
        axis.text.x = element_text(size = 30,family = "Arial"),
        axis.line = element_line(size = 2))

ggsave("SETBP1_epigenomics/pipeline/plots/Corr_ATAC_H3K7ac_NPC_D868D.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 300, height = 250, units = "mm", dpi = 300, limitsize = TRUE)

ATAC_H3K27ac_NPC_D868N <- ATAC_H3K27ac_NPC_D868  %>% 
  dplyr::filter(NPC_D868N_ATAC_peaks==1 & NPC_D868N_H3K27ac_peaks==1) %>% 
  dplyr::select(NPC_D868N,NPC_D868N_H3K27ac)

ggplot(ATAC_H3K27ac_NPC_D868N, aes(x = log2(NPC_D868N), y = log2(NPC_D868N_H3K27ac))) + 
  geom_point() + 
  geom_abline(colour = "blue")+
  stat_cor(method = "pearson", label.x = 7, label.y = 13, size=6)+
  xlab('ATAC log2(RPKM) NPCs D868N') +
  ylab('H3K27ac log2(RPKM) NPCs D868N')+
  theme_classic()+
  ylim(0,15)+
  xlim(0,13)+
  theme(axis.text.y = element_text(size = 30,family = "Arial"),
        axis.title.y = element_text(size = 28,family = "Arial"),
        axis.title.x = element_text(size = 28,family = "Arial"),
        axis.text.x = element_text(size = 30,family = "Arial"),
        axis.line = element_line(size = 2))

ggsave("SETBP1_epigenomics/pipeline/plots/Corr_ATAC_H3K7ac_NPC_D868N.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 300, height = 250, units = "mm", dpi = 300, limitsize = TRUE)


#Fig.S1 g PCA of ATAC samples NPCs, IPSCs, Zebrafish SGS and inducible SET overexpression lines

