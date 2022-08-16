library(tidyverse)
library(IRanges)
library(dplyr)
library(plyranges)
library(diffloop)
library(LSD)


#Fig.3 c plot of relation between interaction frequency and distance in NPCs D868D

#Extracting interaction frequency and distance relationship form HiC maps at 50kb resolution for chr3 

bin_IF_NPC_D868D <- read_delim("~/Setbp1_Gdrive/setbp1/pipeline/annotations/mega (2f06dcec)/NPC-SETBP1_D868D/straw_norm.tsv.gz",
                                                       delim="\t",col_names = T) %>% 
dplyr::mutate(distance=abs(bin2_start-bin1_end))

bin_IF_NPC_D868D_chr3 <- bin_IF_NPC_D868D %>% 
  dplyr::filter(bin1_chr==c("chr3") & bin2_chr==c("chr3"))%>% 
  dplyr::group_by(distance) %>% 
  dplyr::summarise(IF=mean(IF))

bin_IF_NPC_D868D_chr3$condition <- "NPCs D868D"

ggplot(data=bin_IF_NPC_D868D_chr3, aes(x=distance, y=log2(IF), group=1)) +
  geom_line()+ theme_classic()
ggsave("SETBP1_epigenomics/pipeline/plots/distance_IF_D868D.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 145, units = "mm", dpi = 300, limitsize = TRUE)




bin_IF_NPC_D868N <- read_delim("~/Setbp1_Gdrive/setbp1/pipeline/annotations/mega (2f06dcec)/NPC-SETBP1_D868N/straw_norm.tsv.gz",
                               delim="\t",col_names = T) %>% 
  dplyr::mutate(distance=abs(bin2_start-bin1_end))

bin_IF_NPC_D868N_chr3 <- bin_IF_NPC_D868N %>% 
  dplyr::filter(bin1_chr==c("chr3") & bin2_chr==c("chr3"))%>% 
  dplyr::group_by(distance) %>% 
  dplyr::summarise(IF=mean(IF))

bin_IF_NPC_D868N_chr3$condition <- "NPCs D868N"

ggplot(data=bin_IF_NPC_D868N_chr3, aes(x=distance, y=log2(IF), group=1)) +
  geom_line()+ theme_classic()
ggsave("SETBP1_epigenomics/pipeline/plots/distance_IF_D868N.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


#plotting the relationship between interaction frequency and distance

chr3 <- rbind.data.frame(bin_IF_NPC_D868D_chr3,bin_IF_NPC_D868N_chr3) %>% 
  dplyr::filter(distance>=1000000) %>% 
  dplyr::filter(distance<=10000000) %>% 
  dplyr::mutate(distance=distance/1000000)
  

chr3$IF <- rescale(chr3$IF, from = c(0, 7060.37), to = c(0, 1))

ggplot(data=chr3, aes(x=distance, y=log2(IF), group=condition, color=condition)) +
  geom_line() +
  scale_color_manual(values=c("#006d2c","#74c476"))+
  theme_classic() +xlab('Distance (Mbp)') + 
  ylab('Log2 (Interaction Frequency)')+                                      # Change decimal comma / point  
  theme(panel.border = element_rect())+
  xlim(1,10)+
  ylim(-11,-6)+
  theme_classic()+ theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial", colour = "black"),
        axis.text.y = element_text(size = 18,family = "Arial", colour = "black"),
        axis.title.y = element_text(size = 22,family = "Arial"),
        axis.title.x = element_text(size = 22,family = "Arial"),
        axis.line = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))

ggsave("SETBP1_epigenomics/pipeline/plots/distance_IF_NPC_CHR3.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 130, height = 125, units = "mm", dpi = 300, limitsize = TRUE)

#Fig.3 d Scatter plot loop strength and pie charts differential loops


#Calculate NPCs loop strength of NPCs D868D & D868N in NPCs D868D significant loops 

Loop_NPC_D868D_5kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_5000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868D_10kb <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/enriched_pixels_10000.bedpe",
                                 delim="\t",col_names = T)

Loop_NPC_D868D_inD868N_5kb <- read_delim("SETBP1_epigenomics/HiC/NPC_D868D_vs_NPC_D868N//hiccupsdiff/file2/requested_list_5000.bedpe",
                                         delim="\t",col_names = T)

Loop_NPC_D868D_inD868N_10kb <- read_delim("SETBP1_epigenomics/HiC/NPC_D868D_vs_NPC_D868N/hiccupsdiff/file2/requested_list_10000.bedpe",
                                         delim="\t",col_names = T)

Loop_NPC_D868D_merged <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/hiccups/merged_loops.bedpe",
                                    delim="\t",col_names = T)

Loop_NPC_D868D_merged <- Loop_NPC_D868D_merged[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2)

Loop_NPC_D868D_5kb <- Loop_NPC_D868D_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868D=12,
                expectedBL_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868D,
                expectedBL_D868D)

Loop_NPC_D868D_inD868N_5kb <- Loop_NPC_D868D_inD868N_5kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868N=12,
                expectedBL_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868N,
                expectedBL_D868N)

Loop_NPC_D868D_5kb$chr2 <- as.character(Loop_NPC_D868D_5kb$chr2)
Loop_NPC_D868D_inD868N_5kb$chr2 <- as.character(Loop_NPC_D868D_inD868N_5kb$chr2)
Loop_NPC_D868D_vs_D868N_5kb <- inner_join(Loop_NPC_D868D_5kb,Loop_NPC_D868D_inD868N_5kb)


Loop_NPC_D868D_10kb <- Loop_NPC_D868D_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868D=12,
                expectedBL_D868D=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868D,
                expectedBL_D868D)

Loop_NPC_D868D_inD868N_10kb <- Loop_NPC_D868D_inD868N_10kb[-c(1),] %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                chr2=4,
                start2=5,
                end2=6,
                observed_D868N=12,
                expectedBL_D868N=13) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                observed_D868N,
                expectedBL_D868N)

Loop_NPC_D868D_vs_D868N_10kb <- inner_join(Loop_NPC_D868D_10kb,Loop_NPC_D868D_inD868N_10kb)


Loop_NPC_D868D_vs_D868N_5kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_NPC_D868D_vs_D868N_5kb)
Loop_NPC_D868D_vs_D868N_10kb_merged <- inner_join(Loop_NPC_D868D_merged,Loop_NPC_D868D_vs_D868N_10kb)

Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops <- rbind.data.frame(Loop_NPC_D868D_vs_D868N_5kb_merged,Loop_NPC_D868D_vs_D868N_10kb_merged) %>% 
  dplyr::mutate(intensity_D868D=observed_D868D/expectedBL_D868D) %>% 
  dplyr::mutate(intensity_D868N=observed_D868N/expectedBL_D868N) %>% 
  dplyr::mutate(fold=foldchange(intensity_D868N,intensity_D868D)) %>% 
  dplyr::mutate(chr2=chr1)
  

Loop_NPC_D868D_vs_D868N_fc_minus2 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=-2) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_minus2.bedpe")

Loop_NPC_D868D_vs_D868N_fc_minus1.5 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_minus1.5.bedpe")


Loop_NPC_D868D_vs_D868N_fc_plus2 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold>=2)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_minus2.bedpe")

Loop_NPC_D868D_vs_D868N_fc_plus1.5 <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold>=1.5)%>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_plus1.5.bedpe")


Loop_NPC_D868D_vs_D868N_fc_unchanged <- Loop_NPC_D868D_vs_D868N_merged_Ctrl_loops %>% 
  dplyr::filter(fold<=1.5)%>% 
  dplyr::filter(fold>=-1.5) %>% 
  write_tsv("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_vs_D868N_fc_unchanged1.5.bedpe")


#Fig.3 f g h  Heatmap of RNA levels in loop associated genes

#Load RNA-seq data of NPCs

RNA_seq <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/RNA_seq/Results_pipeline_hg38/normalized_counts.tsv") %>% 
  dplyr::select(Gene,Geneid,NPCs_D868D_1,NPCs_D868D_2,NPCs_D868D_4,NPCs_D868D_6,NPCs_D868N_4,NPCs_D868N_5,NPCs_D868N_7, NPCs_D868N_8) %>% 
  dplyr::mutate(NPC_D868D_RNA=(NPCs_D868D_1+NPCs_D868D_2+NPCs_D868D_6+NPCs_D868D_6)/3,
                NPC_D868N_RNA=(NPCs_D868N_4+NPCs_D868N_5+NPCs_D868N_7+NPCs_D868N_8)/3) %>% 
  dplyr::rename(SYMBOL=1)

# Annotate Loop anchors with Gene annotation, find genes in 10KB ranges 

Loop_NPC_D868D_Annotate <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_Annotate",
                                      delim="\t", col_names = T) %>% 
  dplyr::select(chr...1,start...2,end...3,fold...29,distanceToTSS...17,SYMBOL...19,chr...22,start...23,end...24,distanceToTSS...38,SYMBOL...40)
  
Loop_NPC_D868D_Annotate_1 <- Loop_NPC_D868D_Annotate %>% 
  dplyr::select(chr...1,start...2,end...3,fold...29) %>% 
  dplyr::rename(chr=1,start=2,end=3)

Loop_NPC_D868D_Annotate_2 <- Loop_NPC_D868D_Annotate %>% 
  dplyr::select(chr...22,start...23,end...24,fold...29)%>% 
  dplyr::rename(chr=1,start=2,end=3)

Loop_NPC_D868D_Anchors <- rbind.data.frame(Loop_NPC_D868D_Annotate_1,Loop_NPC_D868D_Annotate_2)%>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

TSS <- TSS.human.GRCh38 %>% 
  addchr() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  dplyr::rename(Geneid=1) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


Loop_NPC_D868D_Anchors_TSS <- annotatePeakInBatch(TSS, 
                                                  AnnotationData=Loop_NPC_D868D_Anchors, 
                                                                 output="shortestDistance",
                                                                 PeakLocForDistance="middle",
                                                                 multiple = FALSE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(distance=abs(distancetoFeature))


Loop_NPC_D868D_Anchors_TSS_RNA_seq <- inner_join(Loop_NPC_D868D_Anchors_TSS,RNA_seq)

Loop_NPC_D868D_Anchors_TSS_RNA_seq_genes <- Loop_NPC_D868D_Anchors_TSS_RNA_seq %>% 
  dplyr::select(SYMBOL) %>% 
  unique()

Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff <- inner_join(Loop_NPC_D868D_Anchors_TSS,RNA_seq) %>% 
  dplyr::filter(fold...29 <= -1.5) %>% 
  dplyr::mutate(distance=grangesLoop.end-grangesGene.start) %>% 
  dplyr::mutate(distance=abs(distance))

Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes <- Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff %>% 
  dplyr::select(SYMBOL) %>% 
  unique() %>% 
  write_tsv("SETBP1_epigenomics/HiC/Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes")

#Heatmap of associated genes 

RNA_seq <- read_tsv("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/RNA_seq/Results_pipeline_hg38/normalized_counts.tsv") %>% 
  dplyr::rename(SYMBOL=1)


metadata = data.frame(samples = rep(as.character(c(1, 2, 3, 4)), 4, 8),
                      condition = str_split_fixed(names(RNA_seq)[c(9:10, 12, 14, 19:22)], "_", 3)[,2],
                      row.names = names(RNA_seq)[c(9:10, 12, 14, 19:22)],
                      stringsAsFactors = T)

suppressMessages(library(RColorBrewer))
suppressMessages(library("viridis"))

annotation_column <- metadata[,1:(dim(metadata)[2])]
mycolors_s <- as.vector(polychrome(4)); names(mycolors_s) = levels(annotation_column$samples)
mycolors_c <- c("#006d2c", "#74c476"); names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(samples = mycolors_s, condition=mycolors_c)
crp <- colorRampPalette(c('blue','white','red'))
colors = crp(255)


RNA_seq <- RNA_seq[c(2:3,5,7,12:15)]

Expression_distribution_in_loop_distance <- inner_join(Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes,RNA_seq)

to_H <- Expression_distribution_in_loop_distance[c(2:3,5,7,12:15)]
to_H <- to_H[rowSums(to_H) > 0,]

to_H <- to_H %>%
  dplyr::rename("NPCs D868D 1"=1,
                "NPCs D868D 2"=2,
                "NPCs D868D 3"=3,
                "NPCs D868D 4"=4,
                "NPCs D868N 1"=5,
                "NPCs D868N 2"=6,
                "NPCs D868N 3"=7,
                "NPCs D868N 4"=8)
  

pheatmap(as.matrix(to_H), annotation_col = annotation_column,
         annotation_colors = ann_colors, scale = "row", col=colors, 
         show_rownames = FALSE)

Expression_distribution_in_loop_distance <- inner_join(Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes,RNA_seq)  %>% 
  dplyr::select(NPC_D868D_RNA,NPC_D868N_RNA) %>% 
  dplyr::rename("NPC_D868D"=1,
                "NPC_D868N"=2) %>% 
  gather(key=Group, value=norm_counts, "NPC_D868D","NPC_D868N")
  

ggplot(Expression_distribution_in_loop_distance) +
  aes(x = Group, y = log2(norm_counts), fill = Group, color = Group) +
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
    ref.group = "NPC_D868D",
    label.y = 17  #posizione p value
  )

