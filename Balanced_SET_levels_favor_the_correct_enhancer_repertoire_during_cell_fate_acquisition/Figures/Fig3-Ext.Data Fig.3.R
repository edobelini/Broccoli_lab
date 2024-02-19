library(tidyverse)
library(IRanges)
library(dplyr)
library(plyranges)
library(diffloop)
library(LSD)
library(hictoolsr)
library(plotgardener)
library(pheatmap)

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


  #Gene Ontology on loop genes 

  gostres <- gost(query = unique(Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes$SYMBOL),
                                      organism = "hsapiens",
                                      evcodes = TRUE,
                                      significant = TRUE,
                                      correction_method = "fdr",
                                      user_threshold = 0.05 , sources = c("GO:BP"))

GO <- as.data.frame(gostres$result)
GO$Perc_of_enrichment <- GO$intersection_size / GO$term_size *100
GO$Condition <- "All"

GO <- GO %>% 
  dplyr::mutate(log10_pvalue=-log10(p_value)) %>% 
  dplyr::rename("Percentage of Enrichment"=17,
                "-log10 Pvalue"=19) 

write_tsv(GO,"SETBP1_epigenomics/HiC/GO_Loop_NPC_D868D_Anchors_TSS_RNA_seq_diff_genes")

names(GO)GO$`Percentage of Enrichment`
ggplot(data = GO[c(1:10, 12, 14, 16),], aes(x = Condition, y = reorder(term_name,`Percentage of Enrichment`), color = `-log10 Pvalue`, size = `Percentage of Enrichment`)) +
  geom_point(stroke = 1)+
  scale_color_gradient2(low = "grey", mid = "orange", high = "red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 18, color = 'black',family = "Arial"),
        axis.text.y = element_text(size = 18, color = 'black', hjust = 1, family="Arial"),
        axis.title = element_text(size = 18, color = 'black', hjust = 1, family="Arial"),
        legend.text = element_text(size = 12, color = 'black', hjust = 1, family="Arial"),
        legend.title = element_text(size = 14, color = 'black', hjust = 1, family="Arial"),
        panel.background = element_rect(fill = "white",colour = "grey60", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey40"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
        title = element_text(size = 20, color = 'black', hjust = 1, family="Arial"))+ 
  coord_fixed(ratio = 1)+
  ggtitle("GO Biological Processes")+
  ylab("") +
  xlab("") 

ggsave(filename = "SETBP1_epigenomics/pipeline/plots/GO_Gene_loop_decreased.png", 
       scale = 1, width = 60, height = 115, units = "mm", dpi = 300, limitsize = TRUE)


#Fig.3 Extended Data C Contact domains legth and number 


TAD_NPC_D868D <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868D/mega/arrowhead/5000_blocks.bedpe",
                            delim="\t", col_names = T) %>%
  dplyr::rename(chr=1,
                start=2,
                end=3)

TAD_NPC_D868D <- TAD_NPC_D868D[-c(1),] %>% 
  dplyr::select(chr,start,end) %>%
  dplyr::mutate(end-start) %>% 
  dplyr::rename(TAD_length=4) 

  TAD_NPC_D868D$condition <- "NPCs D868D"

TAD_NPC_D868N <- read_delim("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/NPC_D868N/mega/arrowhead/5000_blocks.bedpe",
                            delim="\t", col_names = T) %>% 
  dplyr::rename(chr=1,
         start=2,
         end=3)

TAD_NPC_D868N <- TAD_NPC_D868N[-c(1),] %>% 
  dplyr::select(chr,start,end) %>%
  dplyr::mutate(end-start) %>% 
  dplyr::rename(TAD_length=4)


TAD_NPC_D868N$condition <- "NPCs D868N"

genotype <- c("NPCs D868D","NPCs D868N")
value <- c( 6343,
            5369)
TADs <- data.frame(genotype,value) 

data1 <- TADs                                                 
data1$genotype <- factor(data1$genotype,                                   
                         levels = c("NPCs D868D","NPCs D868N"))

ggplot(data1, aes(y=value, x=genotype, fill=genotype)) + 
  geom_bar(stat="identity",width = 0.8,
           size = 1,position = position_dodge(width = 0.2))+
  xlab("")+
  ylab("Contact Domain Number")+
  scale_fill_manual(values = c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial", colour = "black"),
        axis.text.y = element_text(size = 18,family = "Arial", colour = "black"),
        axis.title.y = element_text(size = 22,family = "Arial"),
        axis.line = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))

ggsave("SETBP1_epigenomics/pipeline/plots/Compartments_number_NPC.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)

TAD_length_NPC <- rbind.data.frame(TAD_NPC_D868D,TAD_NPC_D868N)


ggplot(TAD_length_NPC) +
  aes(x = condition, y = TAD_length, fill = condition, color = condition) +
  geom_violin() +
  scale_fill_manual(values = c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  ggthemes::theme_base() +
  xlab('') + 
  ylab('Contact Domain Length (bp)')+                                      # Change decimal comma / point  
  scale_y_continuous(labels = scales::comma_format(big.mark = ".",
                                           decimal.mark = ","))+
  theme(panel.border = element_rect())+
  theme_classic()+ theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial", colour = "black"),
        axis.text.y = element_text(size = 18,family = "Arial", colour = "black"),
        axis.title.y = element_text(size = 22,family = "Arial"),
        axis.line = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))+
  stat_summary(fun.data = "mean_sdl", geom = "crossbar",
               mult=1,fill="white",
               colour = "black", width = 0.05) +
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    ref.group = "Neu D868D",
    label.y =   #posizione p value
  )
  
ggsave("SETBP1_epigenomics/pipeline/plots/Compartments_length_NPC.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)

#Fig.3 Extended Data D Loop number 


genotype <- c("NPCs D868D","NPCs D868N")
value <- c( 14168,
            12815)
TADs <- data.frame(genotype,value) 

data1 <- TADs                                                 
data1$genotype <- factor(data1$genotype,                                   
                         levels = c("NPCs D868D","NPCs D868N"))

ggplot(data1, aes(y=value, x=genotype, fill=genotype)) + 
  geom_bar(stat="identity",width = 0.8,
           size = 1,position = position_dodge(width = 0.2))+
  xlab("")+
  ylab("Chromatin Loops Number")+
  scale_fill_manual(values = c("#006d2c","#74c476"))+
  scale_color_manual(values=c("#006d2c","#74c476")) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 18,family = "Arial", colour = "black"),
        axis.text.y = element_text(size = 18,family = "Arial", colour = "black"),
        axis.title.y = element_text(size = 22,family = "Arial"),
        axis.line = element_line(size = 1),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))

ggsave("SETBP1_epigenomics/pipeline/plots/Loops_number_NPC.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 150, height = 145, units = "mm", dpi = 300, limitsize = TRUE)


#Fig.3 Extended Data E APA plots in different loops subset 


#Differential loop APA plotting 

## Create divergent matrix ####
m <- matrix(data = rnorm(n = 21*21, mean = 0, sd = 2), nrow = 21, ncol = 21)
m1 <- read.csv("/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/APA/NPC_D868D_diff/5000/gw/APA.txt", header=FALSE) %>% 
  as.matrix()
m2 <- read.csv("/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/APA/NPC_D868D_diff/D868N/5000/gw/APA.txt", header=FALSE) %>% 
  as.matrix()
log2(m1)
## Define parameters
p <- pgParams(width = 3, height = 3, default.units = "inches")

## Create page
pageCreate(params = p)

## Plot apa
plot <- plotApa(apa =m1,
                x = p$width/2, y = p$height/2,
                width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
                palette = colorRampPalette(c("blue","white","red")), zrange = c(0,100000)) #plot image was cut & paste by R studio

## Annotate legend
annoHeatmapLegend(plot = plot,
                  x = 4, y = 0.75, width = 0.1, height = 0.75)

## Define parameters
p <- pgParams(width = 3, height = 3, default.units = "inches")

## Create page
pageCreate(params = p)


## Plot apa
plot <- plotApa(apa =m2,
                x = p$width/2, y = p$height/2,
                width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
                palette = colorRampPalette(c("blue","white","red")), zrange = c(0,100000)) #plot image was cut & paste by R studio

#Non-Differential loop APA plotting

m <- matrix(data = rnorm(n = 21*21, mean = 0, sd = 2), nrow = 21, ncol = 21)
m1 <- read.csv("/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/APA/NPC_D868D_no_diff/5000/gw/APA.txt", header=FALSE) %>% 
  as.matrix()
m2 <- read.csv("/home/zaghi/Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Hi-C/hg38_HiC/APA/NPC_D868D_no_diff/D868N/5000/gw/APA.txt", header=FALSE) %>% 
  as.matrix()
log2(m1)
## Define parameters
p <- pgParams(width = 3, height = 3, default.units = "inches")

## Create page
pageCreate(params = p)

## Plot apa
plot <- plotApa(apa =m1,
                x = p$width/2, y = p$height/2,
                width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
                palette = colorRampPalette(c("blue","white","red")), zrange = c(0,100000)) #plot image was cut & paste by R studio

## Annotate legend
annoHeatmapLegend(plot = plot,
                  x = 4, y = 0.75, width = 0.1, height = 0.75)

## Define parameters
p <- pgParams(width = 3, height = 3, default.units = "inches")

## Create page
pageCreate(params = p)


## Plot apa
plot <- plotApa(apa =m2,
                x = p$width/2, y = p$height/2,
                width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
                palette = colorRampPalette(c("blue","white","red")), zrange = c(0,100000)) #plot image was cut & paste by R studio


#Fig.3 Extended Data F Loop length comparison (differential loop vs random subset of unchanged loops of the same number)


#load total set of loops with annotation 



Loop_NPC_D868D_Annotate <- read_delim("SETBP1_epigenomics/pipeline/Peaks/Loop_NPC_D868D_Annotate",
                                      delim="\t", col_names = T) %>% 
  dplyr::rename(chr1=1,
                start1=2,
                end1=3,
                intensity_D868D1=6,
                intensity_D868N1=7,
                foldchange1=8,
                distanceToTSS1=17,
                SYMBOL1=19,
                Feature1=21,
                chr2=22,
                start2=23,
                end2=24,
                distanceToTSS2=38,
                SYMBOL2=40,
                Feature2=42) %>% 
  dplyr::select(chr1,
                start1,
                end1,
                intensity_D868D1,
                intensity_D868N1,
                foldchange1,
                distanceToTSS1,
                SYMBOL1,
                Feature1,
                chr2,
                start2,
                end2,
                distanceToTSS2,
                SYMBOL2,
                Feature2)


#create annotated dataset of unchanged loops and decreased loops

Loop_NPC_D868D_Annotate.bedpe <- Loop_NPC_D868D_Annotate %>% 
  dplyr::mutate(lenght=end2-start1)

Loop_NPC_D868D_Annotate.bedpe$lenght <- abs(Loop_NPC_D868D_Annotate.bedpe$lenght)

mean(Loop_NPC_D868D_Annotate.bedpe$lenght)


Loop_NPC_D868D_Annotate.bedpe_decrease <- Loop_NPC_D868D_Annotate.bedpe %>% 
  dplyr::filter(foldchange1 <= -1.5)

Loop_NPC_D868D_Annotate.bedpe_decrease$condition <- "Decreased"

Loop_NPC_D868D_Annotate.bedpe_unchanged <- Loop_NPC_D868D_Annotate.bedpe %>% 
  dplyr::filter(foldchange1 >= -1.5) %>% 
  dplyr::filter(foldchange1 <= 1.5)

Loop_NPC_D868D_Annotate.bedpe_unchanged$condition <- "Unchanged"

#Create random subsamples of loops with the same dimensione of the decreased loop dataset and plot 

Loop_NPC_D868D_Annotate.bedpe_unchanged$id <- 1:nrow(Loop_NPC_D868D_Annotate.bedpe_unchanged)

Random_loop_1 <- Loop_NPC_D868D_Annotate.bedpe_unchanged[Loop_NPC_D868D_Annotate.bedpe_unchanged$id %in%  sample(Loop_NPC_D868D_Annotate.bedpe_unchanged$id, 2978),] %>% 
  dplyr::rename(lenght_random1=16)

Random_loop_2 <- Loop_NPC_D868D_Annotate.bedpe_unchanged[Loop_NPC_D868D_Annotate.bedpe_unchanged$id %in%  sample(Loop_NPC_D868D_Annotate.bedpe_unchanged$id, 2978),] %>% 
  dplyr::rename(lenght_random2=16)

Random_loop_3 <- Loop_NPC_D868D_Annotate.bedpe_unchanged[Loop_NPC_D868D_Annotate.bedpe_unchanged$id %in%  sample(Loop_NPC_D868D_Annotate.bedpe_unchanged$id, 2978),] %>% 
  dplyr::rename(lenght_random3=16)

lenght_loop <- data.frame(Decreased=Loop_NPC_D868D_Annotate.bedpe_decrease$lenght,Random1=Random_loop_1$lenght_random1,
                          Random2=Random_loop_2$lenght_random2,Random3=Random_loop_3$lenght_random3)

Violin_plot_lenght_loop <- lenght_loop %>%
  gather(key=Group, value=length, Decreased,Random1,Random2,Random3)

ggplot(Violin_plot_lenght_loop) +
  aes(x = Group, y = length, fill = Group, color = Group) +
  geom_violin() +
  scale_fill_manual(values=c("#404040","#BABABA","#BABABA","#BABABA"))+
  scale_color_manual(values=c("#404040","#BABABA","#BABABA","#BABABA")) +
  ggthemes::theme_base() +
  xlab('') + 
  ylab('Anchors distance in KB')+
  theme(panel.border = element_rect())+                                    # Change decimal comma / point  
  scale_y_continuous(labels = scales::comma_format(big.mark = ".",
                                                   decimal.mark = ","))+
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
    ref.group = "Decreased",
    label.y = 6500000  #posizione p value
  )

ggsave("SETBP1_epigenomics/pipeline/plots/Loop_dimension_random_unchanged.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 195, units = "mm", dpi = 300, limitsize = TRUE)   




