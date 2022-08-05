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


#Fig.1 g Violin plot ATAC in all NPC D868 peaks

Violin_plot_ATAC_NPCs_all <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(NPC_D868D_peaks==1 | NPC_D868N_peaks==1) %>%
  dplyr::select(NPC_D868D,NPC_D868N) %>% 
  dplyr::rename("NPCs D868D"=1,
                "NPCs D868N"=2) %>% 
  gather(key=Group, value=RPKM, "NPCs D868D","NPCs D868N")

t.test(Violin_plot_ATAC_NPCs_Ctrl)

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

t.test(Violin_plot_ATAC_NPCs_Ctrl)

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


