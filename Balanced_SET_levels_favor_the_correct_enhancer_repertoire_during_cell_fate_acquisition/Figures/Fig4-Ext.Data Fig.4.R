library(dplyr)
library(tidyverse)
library(ggpubr)



#Load Table containing all annotation relative to NPCs and Neurons from SGS patient line D868 as obtained and presented in Fig.1


Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) 

#Fig. 4 differential accessibility in development SGS NPCs to NPCs-derived neurons plot (to understand how differential accessibility was calculated refer to the differential accessibility in development DESEQ2 code)


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
