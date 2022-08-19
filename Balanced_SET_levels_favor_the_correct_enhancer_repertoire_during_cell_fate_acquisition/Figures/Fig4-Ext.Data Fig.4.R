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


#Fig 4 D Loop anchors contact strength dynamic changes during Neural differentation

#Loop strength plot NPCs D868D vs Neu D868D (Data calculated as in the Loop_strength_caculation_Neuro_dev.R script)

png(filename = "SETBP1_epigenomics/pipeline/plots/Loop_strenght_all_Neu_D868D_vs_NPC_D868D.png", width = 25, height = 25, res = 720, units = "cm")
heatscatter(x=Loop_NPC_D868D_vs_Neu_D868D_all_loops$intensity_NPC_D868D,y=Loop_NPC_D868D_vs_Neu_D868D_all_loops$intensity_Neu_D868D, 
            colpal="bl2gr2rd", 
            main="", 
            cor=FALSE,
            xlab = "log2(Obs/ExpBL) NPC D868D", 
            ylab = "log2(Obs/ExpBL) Neu D868D",
            ylim = c(-1,15),
            xlim = c(-1,15),
            #only = "x",
            cexplot = 0.6, 
            nrcol = 100,  
            grid = 100)+
  abline(a = 0, b = 1, lwd=3, lty=2)
dev.off()

Loops_1.5fold_piecharts <- data.frame(Modification= c('Decreased','Increased','Unchanged'),
                                      value=c(15.45,28.24,56.31))

ggplot(Loops_1.5fold_piecharts, aes(x="", y=value, fill=Modification)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()
ggsave("SETBP1_epigenomics/pipeline/plots/Loops_1.5fold_piecharts_Neu_vs_NPC_D868D.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 80, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

#Loop strength plot NPCs D868N vs Neu D868N (Data calculated as in the Loop_strength_caculation_Neuro_dev.R script)

png(filename = "SETBP1_epigenomics/pipeline/plots/Loop_strenght_all_Neu_D868N_vs_NPC_D868N.png", width = 25, height = 25, res = 720, units = "cm")
heatscatter(x=Loop_NPC_D868N_vs_Neu_D868N_all_loops$intensity_NPC_D868N,y=Loop_NPC_D868N_vs_Neu_D868N_all_loops$intensity_Neu_D868N, 
            colpal="bl2gr2rd", 
            main="", 
            cor=FALSE,
            xlab = "log2(Obs/ExpBL) NPC D868N", 
            ylab = "log2(Obs/ExpBL) Neu D868N",
            ylim = c(-1,15),
            xlim = c(-1,15),
            #only = "x",
            cexplot = 0.6, 
            nrcol = 100,  
            grid = 100)+
  abline(a = 0, b = 1, lwd=3, lty=2)
dev.off()

Loops_1.5fold_piecharts <- data.frame(Modification= c('Decreased','Increased','Unchanged'),
                                      value=c(6.80,33.23,59.97))

ggplot(Loops_1.5fold_piecharts, aes(x="", y=value, fill=Modification)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()
ggsave("SETBP1_epigenomics/pipeline/plots/Loops_1.5fold_piecharts_Neu_vs_NPC_D868N.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 80, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


#Fig 4 e Loop anchors differential loops anchors gaining strength in development 

#Extract single anchors from loop anchors left and right

Loop_NPC_D868D_vs_Neu_D868D_fc1.5_plus_anchor1 <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_Neu_D868D_all_loops_fc_plus1.5.bedpe",
                                                     delim="\t", col_names = T) %>% 
  dplyr::select(chr1,start1,end1) %>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)



Loop_NPC_D868N_vs_Neu_D868N_fc1.5_plus_anchor1 <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops_fc_plus1.5.bedpe",
                                                             delim="\t", col_names = T) %>% 
  dplyr::select(chr1,start1,end1) %>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

Loop_NPC_D868D_vs_Neu_D868D_fc1.5_plus_anchor2 <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_Neu_D868D_all_loops_fc_plus1.5.bedpe",
                                                             delim="\t", col_names = T) %>% 
  dplyr::select(chr2,start2,end2) %>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)



Loop_NPC_D868N_vs_Neu_D868N_fc1.5_plus_anchor2 <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868N_vs_Neu_D868N_all_loops_fc_plus1.5.bedpe",
                                                             delim="\t", col_names = T) %>% 
  dplyr::select(chr2,start2,end2) %>% 
  dplyr::rename(chr=1,start=2,end=3) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


#Create a unique loop list for D868D & D868N

Loop_anchor_NPC_D868D_vs_Neu_D868D_fc1.5_plus <- rbind.data.frame(Loop_NPC_D868D_vs_Neu_D868D_fc1.5_plus_anchor1,Loop_NPC_D868D_vs_Neu_D868D_fc1.5_plus_anchor2) %>% 
  makeGRangesFromDataFrame() %>% 
  addchr()%>% 
  as.data.frame() %>% 
  unique()

Loop_anchor_NPC_D868N_vs_Neu_D868N_fc1.5_plus <- rbind.data.frame(Loop_NPC_D868N_vs_Neu_D868N_fc1.5_plus_anchor1,Loop_NPC_D868N_vs_Neu_D868N_fc1.5_plus_anchor2) %>% 
  makeGRangesFromDataFrame() %>% 
  addchr() %>% 
  as.data.frame() %>% 
  unique()

#Define different subset of loops anchors for Venn Diagram 

Loop_anchor_Neu_fc1.5_plus_common <- inner_join(Loop_anchor_NPC_D868D_vs_Neu_D868D_fc1.5_plus,
                                                Loop_anchor_NPC_D868N_vs_Neu_D868N_fc1.5_plus, by=c("seqnames","start","end")) %>% 
  unique() %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Common_loop_anchor_Neu_up_in_dev.bedpe")


Loop_anchor_Neu_fc1.5_plus_D868D_only <- anti_join(Loop_anchor_NPC_D868D_vs_Neu_D868D_fc1.5_plus,
                                                   Loop_anchor_NPC_D868N_vs_Neu_D868N_fc1.5_plus, by=c("seqnames","start","end")) %>% 
  unique() %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/loop_anchor_D868D_only_Neu_up_in_dev.bedpe")

Loop_anchor_Neu_fc1.5_plus_D868N_only <- anti_join(Loop_anchor_NPC_D868N_vs_Neu_D868N_fc1.5_plus,
                                                   Loop_anchor_NPC_D868D_vs_Neu_D868D_fc1.5_plus, by=c("seqnames","start","end")) %>% 
  unique() %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/loop_anchor_D868N_only_Neu_up_in_dev.bedpe")


#Fig 4 f Differential regulation genes associated to peaks upregulated in development 


ATAC_peaks_neurons_D868D_HiC <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_HiC") #Tables calculated Enhancers-promoters contact script  

ATAC_peaks_neurons_D868N_HiC <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_HiC")  #Tables calculated Enhancers-promoters contact script 

Coregulated_genes <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/Common_coregulated_genes_SYMBOL") #Coregulated genes calculated in the Enhancers-promoters contact script 
  
  
ATAC_peaks_neurons_D868D_HiC_dev_count <- inner_join(Coregulated_genes,ATAC_peaks_neurons_D868D_HiC,by=c("SYMBOL")) %>% 
  group_by(SYMBOL) %>%
  dplyr::count()

ATAC_peaks_neurons_D868N_HiC_dev_count <- inner_join(Coregulated_genes,ATAC_peaks_neurons_D868N_HiC,by=c("SYMBOL")) %>% 
  group_by(SYMBOL) %>%
  dplyr::count()

ATAC_peaks_neurons_D868D_HiC_dev_coregulated <- inner_join(Coregulated_genes,ATAC_peaks_neurons_D868D_HiC,by=c("SYMBOL")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

ATAC_peaks_neurons_D868N_HiC_dev_coregulated <- inner_join(Coregulated_genes,ATAC_peaks_neurons_D868N_HiC,by=c("SYMBOL")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

common_coregulated_peaks_count_all <- join_overlap_intersect(ATAC_peaks_neurons_D868D_HiC_dev_coregulated,
                                                             ATAC_peaks_neurons_D868N_HiC_dev_coregulated,suffix = c("D868D", "D868N")) %>% 
  as.data.frame() %>% 
  dplyr::filter(SYMBOLD868D == SYMBOLD868N) %>% 
  dplyr::rename(SYMBOL=6) %>% 
  dplyr::select(SYMBOL,peaksD868N) %>% 
  unique() %>% 
  group_by(SYMBOL)%>%
  dplyr::count()

all_counts_all_regions <- full_join(full_join(ATAC_peaks_neurons_D868D_HiC_dev_count,ATAC_peaks_neurons_D868N_HiC_dev_count, by = "SYMBOL"), common_coregulated_peaks_count_all, by = "SYMBOL")
all_counts_all_regions[is.na(all_counts_all_regions)] <- 0 
names(all_counts_all_regions)[2:4] <- c("D868D_count", "D868N_count", "Common_count")

all_counts_all_regions$Unique_D868D <- all_counts_all_regions$D868D_count - all_counts_all_regions$Common_count
all_counts_all_regions$Unique_D868N <- all_counts_all_regions$D868N_count - all_counts_all_regions$Common_count

#Venn diagram differential regulation 

gene_same_regulation <- all_counts_all_regions %>% 
  dplyr::filter(Unique_D868D==0 & Unique_D868N==0)

gene_mixed_regulation <- all_counts_all_regions %>% 
  dplyr::filter(Common_count>0)

gene_diff_regulation <- all_counts_all_regions %>% 
  dplyr::filter(Common_count==0)

Gene_regulation_up_in_dev <- data.frame(
  group = c("Same","Differential"),
  value = c(72,4155)
)

bp<- ggplot(Gene_regulation_up_in_dev, aes(x="", y=value, fill=group))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3"))+ theme_void()
ggsave("SETBP1_epigenomics/pipeline/plots/piecharts_Coregulated_all_in_dev.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 80, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


#Fig 4 h Functional enrichment in Coregulated genes

ATAC_peaks_neurons_D868D_non_coregulated_SYMBOL <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868D_non_coregulated_SYMBOL") #Neu D868D specific genes calculated in the Enhancers-promoters contact scrip

ATAC_peaks_neurons_D868N_non_coregulated_SYMBOL <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/ATAC_peaks_neurons_D868N_non_coregulated_SYMBOL") #Neu D868N specific genes calculated in the Enhancers-promoters contact script

Coregulated_genes <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/Common_coregulated_genes_SYMBOL") #Coregulated genes calculated in the Enhancers-promoters contact script


GO_list <- list(`Neu D868D`=ATAC_peaks_neurons_D868D_non_coregulated_SYMBOL$SYMBOL,
                `Neu D868N`=ATAC_peaks_neurons_D868N_non_coregulated_SYMBOL$SYMBOL,
                `Neu D868D|Neu D868N`=Coregulated_genes$SYMBOL)
library(gprofiler2)
GO_out <- list()
for (i in 1:3) {
  gostres <- gost(query = GO_list[[i]],
                  organism = "hsapiens",
                  evcodes = TRUE,
                  significant = TRUE,
                  correction_method = "fdr",
                  user_threshold = 0.05 , sources = c("GO:BP"))
  
  result <- as.data.frame(gostres$result)
  GO_out[[i]] <- result
  names(GO_out)[i] <- names(GO_list)[i]
}

for (i in 1:3) {
  GO_out[[i]]$Condition <- names(GO_out)[i]
}

GO_all <- bind_rows(GO_out)
GO_all <- GO_all[GO_all$term_name %like% c("neuron"),]

GO_all$Perc_of_enrichment <- GO_all$intersection_size / GO_all$term_size *100
GO_all$Log10_Pvalue <- -log10(GO_all$p_value) 

GO_all <- GO_all %>% 
  dplyr::rename(Condition=17,
                "Percentage of enrichment"=18,
                "-log10 Pvalue"=19)
GO_all <- GO_all[order(GO_all$ `-log10 Pvalue`),]

write_tsv(GO_all,"SETBP1_epigenomics/pipeline/enhancer_promoter/GO_all_coregulated.tsv")

GO_Common <- GO_all %>% 
  dplyr::filter(`-log10 Pvalue`>=5 & Condition==c("Neu D868D|Neu D868N"))

GO_all_2 <- inner_join(GO_Common,GO_all,by="term_name") %>% 
  dplyr::rename(Condition=35,
                "Percentage of enrichment"=36,
                "-log10 Pvalue"=37)

ggplot(data = GO_all_2, aes(x = Condition, y = reorder(term_name, `-log10 Pvalue`), color = `-log10 Pvalue`, size = `Percentage of enrichment`)) +
  geom_point(stroke = 1)+
  ggtitle("GO Biological Processes")+
  theme(axis.title = element_text(size = 35, color = 'black', hjust = 1, family="Arial"))+
  scale_color_gradient2(low = "grey", mid = "orange", high = "red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 18, color = 'black',family = "Arial"),
        axis.text.y = element_text(size = 18, color = 'black', hjust = 1, family="Arial"),
        legend.text = element_text(size = 18, color = 'black', hjust = 1, family="Arial"),
        legend.title = element_text(size = 20, color = 'black', hjust = 1, family="Arial"),
        panel.background = element_rect(fill = "white",colour = "grey60", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey40"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"))+ 
  coord_fixed(ratio = 1)+
  ylab("") +
  xlab("") 
ggplot2::ggsave(filename = "SETBP1_epigenomics/pipeline/plots/GO.png",plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 400, height = 395, units = "mm", dpi = 300, limitsize = TRUE)


#Fig 4 i Heatmap common coregulated genes 

RNA_seq <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/RNA_seq/Results_pipeline_hg38/normalized_counts.tsv") %>% # Load normalized counts 
  dplyr::rename(SYMBOL=1) %>% 
  na.omit()

RNA_seq <- RNA_seq[c(1,3:8)] # Select columns containing only relevant data

Coregulated_genes <- read_tsv("SETBP1_epigenomics/pipeline/enhancer_promoter/Common_coregulated_genes_SYMBOL") #Coregulated genes calculated in the Enhancers-promoters contact script

# Prepare and plot data 

metadata = data.frame(samples = rep(as.character(c(1, 2, 3)), 3, 6),
                      condition = str_split_fixed(names(RNA_seq)[c(2:7)], "_", 3)[,2],
                      row.names = names(RNA_seq)[c(2:7)],
                      stringsAsFactors = T)

suppressMessages(library(RColorBrewer))
suppressMessages(library("viridis"))

annotation_column <- metadata[,1:(dim(metadata)[2])]
mycolors_s <- c("white", "white","white"); names(mycolors_s) = levels(annotation_column$samples)
mycolors_c <- c("#08519c", "#6baed6"); names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(samples = mycolors_s, condition=mycolors_c)
crp <- colorRampPalette(c('blue','white','red'))
colors = crp(255)

Expression_common_coregulated <- inner_join(Coregulated_genes,RNA_seq)

to_H <- Expression_common_coregulated 
to_H <- to_H[2:7]
to_H <- to_H[rowSums(to_H) > 0,]

to_H <- to_H 

h1 <- pheatmap(as.matrix(to_H), annotation_col = annotation_column,
         annotation_colors = ann_colors, scale = "row", col=colors, cluster_cols = T,
         show_rownames = FALSE,annotation_names_col = F,annotation_legend = F,
         labels_col = c("Neu D868D 1","Neu D868D 2", "Neu D868D 3",
                        "Neu D868N 1","Neu D868N 2", "Neu D868N 3"))
h1
png("/home/zaghi/SETBP1_epigenomics/pipeline/plots/heatmap_coregulated.png",pointsize = 1,res=1200,height = 65,width = 35,
    units = "cm")
h1
dev.off()

#Fig. Extended Data 4 a Violin plot and Venn Diagram peaks distribution in clusters 

Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) # Upload tables with peaks annotation and counts (multiBigWigSummary_annotation.sk & multiBigWigSummary_editing.R)

#Violin plot neurons in all Peaks NPCs derived neurons

Violin_plot_ATAC_Neus_all <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(Neu_D868D_peaks==1 | Neu_D868N_peaks==1) %>%
  dplyr::select(Neu_D868D,Neu_D868N)  %>% 
  dplyr::rename("Neu D868D"=1,
                "Neu D868N"=2) %>% 
  gather(key=Group, value=RPKM, "Neu D868D","Neu D868N")

t.test(Violin_plot_ATAC_Neus_all)

ggplot(Violin_plot_ATAC_Neus_all) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#08519c","#6baed6"))+
  scale_color_manual(values=c("#08519c","#6baed6")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial")) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    ref.group = "Neu D868D",
    label.y = 18  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_Neu_D868_all_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE) 

#Violin plot neurons in D868D Peaks NPCs derived neurons

Violin_plot_ATAC_Neu_Ctrl <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(Neu_D868D_peaks==1) %>%
  dplyr::select(Neu_D868D,Neu_D868N) %>% 
  dplyr::rename("Neu D868D"=1,
                "Neu D868N"=2) %>% 
  gather(key=Group, value=RPKM, "Neu D868D","Neu D868N") 

t.test(Violin_plot_ATAC_Neu_Ctrl)

ggplot(Violin_plot_ATAC_Neu_Ctrl) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#08519c","#6baed6"))+
  scale_color_manual(values=c("#08519c","#6baed6")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")+
  theme_classic()+ theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial")) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    ref.group = "Neu D868D",
    label.y = 18  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_Neu_D868D_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)


# Venn Diagram clusters in NPCs derived Neu D868D 


Neu_D868D_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Neu_D868D_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster1.bed")

Neu_D868D_cluster1_annotate <- ChIPseeker::annotatePeak(Neu_D868D_cluster1,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster1_annotate.tsv")


table(Neu_D868D_cluster1_annotate$Feature)

Neu_D868D_cluster1_pie <- data.frame(table(Neu_D868D_cluster1_annotate$Feature))

Neu_D868D_cluster1_pie$percentage <- prop.table(Neu_D868D_cluster1_pie$Freq)*100

bp<- ggplot(Neu_D868D_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



Neu_D868D_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_2"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Neu_D868D_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster2.bed")

Neu_D868D_cluster2_annotate <- ChIPseeker::annotatePeak(Neu_D868D_cluster2,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster2_annotate.tsv")


Neu_D868D_cluster2_pie <- data.frame(table(Neu_D868D_cluster2_annotate$Feature))

Neu_D868D_cluster2_pie$percentage <- prop.table(Neu_D868D_cluster2_pie$Freq)*100

bp<- ggplot(Neu_D868D_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


Neu_D868D_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_3"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Neu_D868D_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster3.bed")



Neu_D868D_cluster3_annotate <- ChIPseeker::annotatePeak(Neu_D868D_cluster3,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster3_annotate.tsv")

Neu_D868D_cluster3_annotate_SYMBOL <- read_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster3_annotate.tsv") %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster3_annotate_SYMBOL.tsv")


Neu_D868D_cluster3_pie <- data.frame(table(Neu_D868D_cluster3_annotate$Feature))

Neu_D868D_cluster3_pie$percentage <- prop.table(Neu_D868D_cluster3_pie$Freq)*100

bp<- ggplot(Neu_D868D_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/Neu_D868D_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

#Fig. Extended Data 4 b Violin plot and Venn Diagram peaks distribution in clusters 

Open_chromatin_all_ATAC_SGS <- read_tsv("SETBP1_epigenomics/pipeline/Peaks/multiBigwigSummary_ATAC_SETBP1_OpenChromatin_table_annotate", col_names = T) # Upload tables with peaks annotation and counts (multiBigWigSummary_annotation.sk & multiBigWigSummary_editing.R)

#Violin plot neurons in all Peaks IPSCs derived neurons

Violin_plot_ATAC_dir_Neus_all <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(dir_Neu_D868D_peaks==1 | dir_Neu_D868N_peaks==1) %>%
  dplyr::select(dir_Neu_D868D,dir_Neu_D868N)  %>% 
  dplyr::rename("dir Neu D868D"=1,
                "dir Neu D868N"=2) %>% 
  gather(key=Group, value=RPKM, "dir Neu D868D","dir Neu D868N")

t.test(Violin_plot_ATAC_dir_Neus_all)

ggplot(Violin_plot_ATAC_dir_Neus_all) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#08519c","#6baed6"))+
  scale_color_manual(values=c("#08519c","#6baed6")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial")) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    ref.group = "dir Neu D868D",
    label.y = 18  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_Neu_D868_all_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE) 

#Violin plot neurons in D868D Peaks IPSCs derived neurons

Violin_plot_ATAC_dir_Neu_Ctrl <- Open_chromatin_all_ATAC_SGS %>%
  dplyr::filter(dir_Neu_D868D_peaks==1) %>%
  dplyr::select(dir_Neu_D868D,dir_Neu_D868N) %>% 
  dplyr::rename("dir Neu D868D"=1,
                "dir Neu D868N"=2) %>% 
  gather(key=Group, value=RPKM, "dir Neu D868D","dir Neu D868N") 

t.test(Violin_plot_ATAC_dir_Neu_Ctrl)

ggplot(Violin_plot_ATAC_dir_Neu_Ctrl) +
  aes(x = Group, y = log2(RPKM), fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#08519c","#6baed6"))+
  scale_color_manual(values=c("#08519c","#6baed6")) +
  ggthemes::theme_base() +
  xlab('') + 
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")+
  theme_classic()+ theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 18,family = "Arial")) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05)+
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    ref.group = "Neu D868D",
    label.y = 18  #posizione p value
  )
ggsave("SETBP1_epigenomics/pipeline/plots/Violin_plot_ATAC_dir_Neu_D868D_peaks.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 120, height = 115, units = "mm", dpi = 300, limitsize = TRUE)


# Venn Diagram clusters in IPSCs derived Neu D868D 

dir_Neu_D868D_cluster1 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_dir_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_1"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Neu_D868D_cluster1,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster1.bed")

dir_Neu_D868D_cluster1_annotate <- ChIPseeker::annotatePeak(dir_Neu_D868D_cluster1,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster1_annotate.tsv")


table(dir_Neu_D868D_cluster1_annotate$Feature)

dir_Neu_D868D_cluster1_pie <- data.frame(table(dir_Neu_D868D_cluster1_annotate$Feature))

dir_Neu_D868D_cluster1_pie$percentage <- prop.table(dir_Neu_D868D_cluster1_pie$Freq)*100

bp<- ggplot(dir_Neu_D868D_cluster1_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster1_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)



dir_Neu_D868D_cluster2 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_dir_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_2"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(dir_Neu_D868D_cluster2,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster2.bed")

dir_Neu_D868D_cluster2_annotate <- ChIPseeker::annotatePeak(dir_Neu_D868D_cluster2,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster2_annotate.tsv")


dir_Neu_D868D_cluster2_pie <- data.frame(table(dir_Neu_D868D_cluster2_annotate$Feature))

dir_Neu_D868D_cluster2_pie$percentage <- prop.table(dir_Neu_D868D_cluster2_pie$Freq)*100

bp<- ggplot(dir_Neu_D868D_cluster2_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster2_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)


dir_Neu_D868D_cluster3 <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/ATAC/pipeline/Compute_Matrix/Peak_centered/Regions/ATAC_in_dir_NeuD868D_ATAC_merge_50_median_Compute_Matrix_3_heatmap.bed")%>% 
  dplyr::filter(deepTools_group==c("cluster_3"))%>% 
  dplyr::rename(chr=1, start=2, end=3) %>% 
  makeGRangesFromDataFrame()

write_bed(Neu_D868D_cluster3,"Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster3.bed")



dir_Neu_D868D_cluster3_annotate <- ChIPseeker::annotatePeak(dir_Neu_D868D_cluster3,
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
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster3_annotate.tsv")

dir_Neu_D868D_cluster3_annotate_SYMBOL <- read_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster3_annotate.tsv") %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster3_annotate_SYMBOL.tsv")


dir_Neu_D868D_cluster3_pie <- data.frame(table(dir_Neu_D868D_cluster3_annotate$Feature))

dir_Neu_D868D_cluster3_pie$percentage <- prop.table(dir_Neu_D868D_cluster3_pie$Freq)*100

bp<- ggplot(dir_Neu_D868D_cluster3_pie, aes(x="", y=Freq, fill=Var1))+
  geom_bar( stat="identity",width=0.5, color="white")

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+
  scale_color_manual(values=c("orange","salmon2","aquamarine3","cyan2","azure3","darkorchid1"))+ theme_void()

ggsave("Setbp1_Gdrive/zaghi_upload/setbp1/Regions/Heatmap_clusters/dir_Neu_D868D_cluster3_pie.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 100, height = 75, units = "mm", dpi = 300, limitsize = TRUE)

