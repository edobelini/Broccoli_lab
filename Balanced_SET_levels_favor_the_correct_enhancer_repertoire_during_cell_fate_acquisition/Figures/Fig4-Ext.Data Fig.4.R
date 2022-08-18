library(dplyr)
library(tidyverse)
library(ggpubr)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(pheatmap)


#Load Table containing all annotation relative to NPCs and Neurons from SGS patient line D868 as obtained and presented in Fig.1


Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) 

#Fig. 4 a differential accessibility in development SGS NPCs to NPCs-derived neurons plot (to understand how differential accessibility was calculated refer to the differential accessibility in development DESEQ2 code)


#NPCs vs Neu D868D

Corrplot_Neu_D868D_dev <- Open_chromatin_all_ATAC_SGS %>% 
  dplyr::filter(Neu_D868D_peaks==1 | NPC_D868D_peaks==1) %>% 
  dplyr::mutate(XY_Neu=Neu_D868D_up_dev-Neu_D868D_down_dev)

ggplot(Corrplot_Neu_D868D_dev) +
  aes(x = log2(NPC_D868D), y = log2(Neu_D868D), colour = as.factor(XY_Neu)) +
  geom_point(size = 1L) +
  ggthemes::theme_base() + 
  theme(legend.position = "none") + 
  geom_abline(slope = 1) +
  scale_color_manual(values=c("#006d2c","grey","#08519c")) +
  labs(x="NPC D868D log2(RPKM)", y="Neu D868D log2(RPKM)")+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"))+
  theme(axis.title.x = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"))
ggsave("SETBP1_Epigenomics/ATAC/Rplots/Corrplot_Neu_D868D_dev_open_Chromatin.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 95, height = 90, units = "mm", dpi = 300, limitsize = TRUE)

#NPCs vs Neu D868N


Corrplot_Neu_D868N_dev <- Open_chromatin_all_ATAC_SGS %>% 
  dplyr::filter(Neu_D868N_peaks==1 | NPC_D868N_peaks==1) %>% 
  dplyr::mutate(XY_Neu=Neu_D868N_up_dev-Neu_D868N_down_dev)

ggplot(Corrplot_Neu_D868N_dev) +
  aes(x = log2(NPC_D868N), y = log2(Neu_D868N), colour = as.factor(XY_Neu)) +
  geom_point(size = 1L) +
  ggthemes::theme_base() + 
  theme(legend.position = "none") + 
  geom_abline(slope = 1) +
  scale_color_manual(values=c("#74c476","grey","#6baed6")) +
  labs(x="NPC D868N log2(RPKM)", y="Neu D868N log2(RPKM)")+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"))+
  theme(axis.title.x = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial"))
ggsave("SETBP1_Epigenomics/ATAC/Rplots/Corrplot_Neu_D868N_dev_open_Chromatin.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 95, height = 90, units = "mm", dpi = 300, limitsize = TRUE)


#Fig. 4 b Violin plot for each different sub group of peaks present in venn diagram in the figure 


#Common Peaks violin plot 


Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(Neu_D868D_up_dev==1 & Neu_D868N_up_dev==1) %>% 
  dplyr::select(NPC_D868D,NPC_D868N,Neu_D868D,Neu_D868N) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "Neu 
 D868D"=3,
                "Neu 
 D868N"=4) %>% 
  gather(key=Group, value=RPKM, "NPCs 
D868D", "NPCs 
D868N", "Neu 
 D868D", "Neu 
 D868N")

my_comparisons <- list(c("NPCs 
D868D","Neu 
 D868D"),c("NPCs 
D868N","Neu 
 D868N"),
                       c("Neu 
 D868D","Neu 
 D868N"))


Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group <- factor(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group, levels = c("NPCs 
D868D","Neu 
 D868D", 
                                                                                                                                                  "NPCs 
D868N","Neu 
 D868N"))


ggplot(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#08519c","#74c476","#6baed6"))+
  scale_color_manual(values=c("#006d2c","#08519c","#74c476","#6baed6")) +
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
    comparisons = my_comparisons,
    label.y = c(17, 18, 20,22)  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_peaks_ATAC_Common_Peaks_neural_development.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


#Up regulated only in Neu D868D


Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(Neu_D868D_up_dev==1 & Neu_D868N_up_dev==0) %>% 
  dplyr::select(NPC_D868D,NPC_D868N,Neu_D868D,Neu_D868N) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "Neu 
 D868D"=3,
                "Neu 
 D868N"=4) %>% 
  gather(key=Group, value=RPKM, "NPCs 
D868D", "NPCs 
D868N", "Neu 
 D868D", "Neu 
 D868N")

my_comparisons <- list(c("NPCs 
D868D","Neu 
 D868D"),c("NPCs 
D868N","Neu 
 D868N"),
                       c("Neu 
 D868D","Neu 
 D868N"))



Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group <- factor(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group, levels = c("NPCs 
D868D","Neu 
 D868D", 
                                                                                                                                                  "NPCs 
D868N","Neu 
 D868N"))


ggplot(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#08519c","#74c476","#6baed6"))+
  scale_color_manual(values=c("#006d2c","#08519c","#74c476","#6baed6")) +
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
    comparisons = my_comparisons,
    label.y = c(12,14, 16,17)  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_peaks_ATAC_D868D_Peaks_neural_development.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)



#Up regulated only in Neu D868N only 

Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(Neu_D868D_up_dev==0 & Neu_D868N_up_dev==1) %>% 
  dplyr::select(NPC_D868D,NPC_D868N,Neu_D868D,Neu_D868N) %>% 
  dplyr::rename("NPCs 
D868D"=1,
                "NPCs 
D868N"=2,
                "Neu 
 D868D"=3,
                "Neu 
 D868N"=4) %>% 
  gather(key=Group, value=RPKM, "NPCs 
D868D", "NPCs 
D868N", "Neu 
 D868D", "Neu 
 D868N")

my_comparisons <- list(c("NPCs 
D868D","Neu 
 D868D"),c("NPCs 
D868N","Neu 
 D868N"),
                       c("Neu 
 D868D","Neu 
 D868N"))



Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group <- factor(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC$Group, levels = c("NPCs 
D868D","Neu 
 D868D", 
                                                                                                                                                  "NPCs 
D868N","Neu 
 D868N"))


ggplot(Box_plot_peaks_ATAC_Common_Peaks_neural_development_HiC) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#006d2c","#08519c","#74c476","#6baed6"))+
  scale_color_manual(values=c("#006d2c","#08519c","#74c476","#6baed6")) +
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
    comparisons = my_comparisons,
    label.y = c(12,14, 16,17)  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_peaks_ATAC_D868N_Peaks_neural_development.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


#Fig 4 C ChromVar Tsne plots (Using ChromVar data calculated in the Relative R script)

#tSNE plots of samples
tsne_plots <- plotDeviationsTsne(dev_Neu, tsne_results, #dev_Neu calculated using the ChromVar script in the Epigenomics folder 
                                 sample_column = "celltype", shiny = FALSE)


tsne_plots

ggsave("SETBP1_epigenomics/pipeline/plots/Tsne_chromVAR_updev.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 115, units = "mm", dpi = 300, limitsize = TRUE)

#tSNE plots of samples with associated deviations score of a specifc TFBS

#NR2F1

tsne_plots_NR2F1 <- plotDeviationsTsne(dev_Neu, tsne_results, annotation = "NR2F1", #dev_Neu calculated using the ChromVar script in the Epigenomics folder
                                 sample_column = "celltype", shiny = FALSE)

tsne_plots_NR2F1

ggsave("SETBP1_epigenomics/pipeline/plots/Tsne_chromVAR_updev_NR2F1.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 115, units = "mm", dpi = 300, limitsize = TRUE)

#EMX1


tsne_plots_EMX1 <- plotDeviationsTsne(dev_Neu, tsne_results, annotation = "EMX1", #dev_Neu calculated using the ChromVar script in the Epigenomics folder
                                 sample_column = "celltype", shiny = FALSE)

tsne_plots_EMX1

ggsave("SETBP1_epigenomics/pipeline/plots/Tsne_chromVAR_updev_NR2F1.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 115, units = "mm", dpi = 300, limitsize = TRUE)





