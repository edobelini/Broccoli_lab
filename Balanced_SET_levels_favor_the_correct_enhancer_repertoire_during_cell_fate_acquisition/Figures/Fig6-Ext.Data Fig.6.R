library(Seurat)
library(dplyr)
library(ArchrR)
library(stringr)
library(ggplot2)
library(openxlsx)
library(pals)
library(RColorBrewer)
library(gprofiler2)
library(ComplexHeatmap)

#Fig.6 e UMAP scMultiome RNA & ATAC


scRNA <- readRDS("Documents/scRNA_Multiome/scRNA_multiome.rds") #load scRNA object 


colors <- kelly(11)
Idents(Multihome)
colors[1] <- "brown" 
p1 = DimPlot(Multihome, reduction = "umap", label = F, pt.size = 0.3, cols = colors,group.by="Label_cluster") + 
  theme_void()+
  theme(legend.text = element_text(color = "black", size  12)) #plot RNA UMAP
ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_RNA.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)


#Add scATAC umap coordinates to Seurat object to plot A

scATAC_Multiome <- loadArchRProject("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/") # Load ArchR project containing scATAC data

barcodes_scATAC <- as.data.frame(scATAC_Multiome@embeddings@listData[["UMAP"]]@listData[["df"]]) %>% 
  rownames_to_column() %>% 
  rename(barcode_scATAC=1)  #extract Barcodes associated with each cells and the UMAP coordinates 

barcodes_scATAC$sample_number <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "-", n = 2)[,2]
barcodes_scATAC$sample <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "#", n = 2)[,1]
barcodes_scATAC$barcode <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "#", n = 2)[,2]
barcodes_scATAC$barcode <- str_split_fixed(barcodes_scATAC$barcode, pattern = "-", n = 2)[,1]

barcodes_scATAC <- barcodes_scATAC %>% 
  rename(sample_numer_ATAC=2)

barcodes_scATAC$sample[barcodes_scATAC$sample == "embryo_mut"] <- "1"
barcodes_scATAC$sample[barcodes_scATAC$sample == "embryo_ctrl"] <- "2"
barcodes_scATAC$sample[barcodes_scATAC$sample == "P2_ctrl"]<- "3"
barcodes_scATAC$sample[barcodes_scATAC$sample == "P2_mut"]<- "4"

barcodes_scATAC$barcode <- paste(barcodes_scATAC$barcode,sep = "-",barcodes_scATAC$sample)
rownames(barcodes_scATAC) <- barcodes_scATAC[,6]
barcodes_scATAC$barcode <- NULL

barcodes_scATAC <- barcodes_scATAC[c(2:3)]

barcodes_scATAC <- barcodes_scATAC %>% 
  dplyr::rename(UMAP_1=1,
                UMAP_2=2)

barcodes_scATAC <- barcodes_scATAC[order(rownames(barcodes)),]
barcodes_scATAC <- barcodes_scATAC[rownames(scRNA@meta.data),] #order the Barcode of scATAC with the same order of scRNA dataset 

scRNA[["cloupe"]] <- CreateDimReducObject(embeddings = scRNA[["umap"]]@cell.embeddings, assay = DefaultAssay(scRNA),
                                          key = "cloupe_") #Insert scATAC UMAP coordinates in Seurat object 

scRNA[["cloupe"]]@cell.embeddings[,1] <- barcodes_scATAC[,1]
scRNA[["cloupe"]]@cell.embeddings[,2] <- barcodes_scATAC[,2]

p1 = DimPlot(Multihome, reduction = "umap", label = F, pt.size = 0.3, cols = colors,group.by="Label_cluster") + 
  theme_void()+
  theme(legend.text = element_text(color = "black", size  12)) #plot ATAC UMAP
ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_ATAC.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)

#Fig.6 g Chromatin accessibility dynamics in pseudotime of neural differentiation 

# Calculate pseudotime trajectory of neural differentiation in control and mutant 

scATAC_Multiome <- loadArchRProject("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR/")

idxPass <- which(scATAC_Multiome$Sample %in% c("embryo_ctrl","P2_ctrl"))
cellsPass <- scATAC_Multiome$cellNames[idxPass]
scATAC_Multiome_ctrl <- scATAC_Multiome[cellsPass, ] #create multiome object with ctrl cells only 

idxPass <- which(scATAC_Multiome$Sample %in% c("embryo_mut","P2_mut"))
cellsPass <- scATAC_Multiome$cellNames[idxPass]
scATAC_Multiome_mut <- scATAC_Multiome[cellsPass, ] #create multiome object with mut cells only 

#Load all the peaks associated with the differentiation trajectory in one list

peaks1 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/AP_RGC") 

peaks2 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/INP") 

peaks3 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/ExN_DL") 

peaks4 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/AP_RGC") 

peaks5 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/INP") 

peaks6 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/ExN_DL") 

peak <- rbind.data.frame(peaks1,peaks2,peaks3,
                         peaks4,peaks5,peaks6) %>% 
  makeGRangesFromDataFrame() %>% 
  IRanges::reduce() #create a unique peak list 

scATAC_Multiome_ctrl <- addPeakSet(
  ArchRProj =  scATAC_Multiome_ctrl,
  peakSet = peak,
  genomeAnnotation = getGenomeAnnotation(scATAC_Multiome),
  force = TRUE
)

scATAC_Multiome_ctrl <- addPeakMatrix(scATAC_Multiome_ctrl, force=T) #Add peakset to ArchR object 

scATAC_Multiome_mut <- addPeakSet(
  ArchRProj =  scATAC_Multiome_mut,
  peakSet = peak,
  genomeAnnotation = getGenomeAnnotation(scATAC_Multiome),
  force = TRUE
)

scATAC_Multiome_mut <- addPeakMatrix(scATAC_Multiome_mut, force=T) #Add peakset to ArchR object 


trajectory_neu <- c("AP_RGC", "INP", "ExN_DL")

scATAC_Multiome_ctrl <- addTrajectory( #addtrajectory calculation to ctrl cells
  ArchRProj = scATAC_Multiome_ctrl, 
  name = "Neurons", 
  groupBy = "Clusters",
  trajectory = trajectory_neu, 
  embedding = "UMAP_2", 
  force = TRUE,
  reducedDims =NULL
)


scATAC_Multiome_mut <- addTrajectory( #addtrajectory calculation to mut cells
  ArchRProj = scATAC_Multiome_mut, 
  name = "Neurons", 
  groupBy = "Clusters",
  trajectory = trajectory_neu, 
  embedding = "UMAP_2", 
  force = TRUE,
  reducedDims =NULL
)

trajPM_ctrl  <- getTrajectory(ArchRProj = scATAC_Multiome_ctrl, name = "Neurons", useMatrix = "PeakMatrix", log2Norm = TRUE) #calculate peaks open during diff. in ctrl.

trajPM_mut <- getTrajectory(ArchRProj = scATAC_Multiome_mut, name = "Neurons", useMatrix = "PeakMatrix", log2Norm = TRUE, ) #calculate peaks open during diff. in mut.

p_trajPM_mut <-  as.data.frame(p_trajPM_mut) %>% 
  rownames_to_column()

p_trajPM_ctrl <-  as.data.frame(p_trajPM_ctrl)%>% 
  rownames_to_column()

colors <- paletteContinuous(set = "solarExtra")

metadata = data.frame(pseudotime = 1:100,
                      stringsAsFactors = T, 
                      row.names=names(p_trajPM_ctrl))

suppressMessages(library(RColorBrewer))
suppressMessages(library("viridis"))

mycolors_s <- paletteContinuous(set = "horizonExtra", n = 100); names(mycolors_s) = unique(metadata$pseudotime)
ann_colors = list(pseudotime=mycolors_s)

p1 <- pheatmap::pheatmap(p_trajPM_ctrl,cluster_rows = F,cluster_cols = F,
                         col=colors,show_rownames = F,
                         show_colnames = F,
                         annotation_col = metadata,
                         annotation_colors = ann_colors) #plot heatmap of peaks in control during differentiation

p2 <- pheatmap::pheatmap(p_trajPM_mut,cluster_rows = F,cluster_cols = F,
                         col=colors,show_rownames = F,
                         show_colnames = F,
                         annotation_col = metadata,
                         annotation_colors = ann_colors) #plot heatmap of peaks in mutant during differentiation

png("SETBP1_epigenomics/pipeline/plots/p_trajPM_ctrl.png",pointsize = 1,res=1200,height = 20,width = 10,
    units = "cm")
p1
dev.off()

png("SETBP1_epigenomics/pipeline/plots/p_trajPM_mut.png",pointsize = 1,res=1200,height = 20,width = 10,
    units = "cm")
p2
dev.off()


p_trajPM_mut_order <- p_trajPM_mut[rownames(p_trajPM_ctrl),] #order peaks of mut in the same order of ctrl.


 p1 <- pheatmap::pheatmap(p_trajPM_mut_order ,cluster_rows = F,cluster_cols = F,
                         col=colors,show_rownames = F,
                         show_colnames = F,
                         annotation_col = metadata,
                         annotation_colors = ann_colors)

png("SETBP1_epigenomics/pipeline/plots/p_trajPM_mut_order.png",pointsize = 1,res=1200,height = 20,width = 10,
    units = "cm")
p1
dev.off()


#Fig.6 h Expression plots during pseudotime of Nrxn1 

p2 <- plotTrajectory(scATAC_Multiome_ctrl, trajectory = "Neurons", colorBy = "GeneScoreMatrix", name = "Nrxn1", continuousSet = "blueYellow")

colors <- paletteContinuous(set = "horizonExtra")
legend_title <- "Pseudotime"

ggplot(p2[[2]][["data"]], aes(x=p2[[2]][["data"]]$x, y=p2[[2]][["data"]]$y, color=p2[[2]][["data"]]$color))+
  geom_point(size = 3) +
  scale_color_gradientn(colours = colors)+
  geom_smooth(color = "black")+
  theme_classic()+
  xlab("PseudoTime")+
  ylab("Nrxn1")+
  ylim(0,25)+
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 20,family = "Arial"),
        axis.title.x = element_text(size = 20,family = "Arial"))
ggsave("SETBP1_epigenomics/pipeline/plots/trajectory_neurons_DL_ctrl_Nrxn1_RNA.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 180, height = 145, units = "mm", dpi = 300, limitsize = TRUE)

p2 <- plotTrajectory(scATAC_Multiome_mut, trajectory = "Neurons", colorBy = "GeneScoreMatrix", name = "Nrxn1", continuousSet = "blueYellow")

colors <- paletteContinuous(set = "horizonExtra")
legend_title <- "Pseudotime"

ggplot(p2[[2]][["data"]], aes(x=p2[[2]][["data"]]$x, y=p2[[2]][["data"]]$y, color=p2[[2]][["data"]]$color))+
  geom_point(size = 3) +
  scale_color_gradientn(colours = colors)+
  geom_smooth(color = "black")+
  theme_classic()+
  xlab("PseudoTime")+
  ylab("Nrxn1")+
  ylim(0,25)+
  theme(axis.text.x = element_text(size = 18,family = "Arial"),
        axis.text.y = element_text(size = 18,family = "Arial"),
        axis.title.y = element_text(size = 20,family = "Arial"),
        axis.title.x = element_text(size = 20,family = "Arial"))
ggsave("SETBP1_epigenomics/pipeline/plots/trajectory_neurons_DL_mut_Nrxn1_RNA.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 180, height = 145, units = "mm", dpi = 300, limitsize = TRUE)



#Fig.6 i Differential enhancers usage plot with degs

#Differential enhancer usage 

#Load peaks associated to each genotype/cluster

files <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/ctrl/", full.names = T)
WT_peaks <- lapply(files, read.table,  header = T)
names(WT_peaks) <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/ctrl/", full.names = F)

files <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/mut/", full.names = T)
MUT_peaks <- lapply(files, read.table,  header = T)
names(MUT_peaks) <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/mut/", full.names = F)


for (i in 1:length(WT_peaks)) {
  WT_peaks[[i]] <- WT_peaks[[i]] %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  MUT_peaks[[i]] <- MUT_peaks[[i]] %>% makeGRangesFromDataFrame(keep.extra.columns = T)
}

common_peaks <- list()

for (i in 1:length(WT_peaks)) {
  common_peaks[[i]] <- subsetByOverlaps(WT_peaks[[i]], MUT_peaks[[i]])
  names(common_peaks)[i] <- names(WT_peaks)[i]
}

for (i in 1:length(WT_peaks)) {
  common_peaks[[i]] <- common_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>%
    as.data.frame()
  
  WT_peaks[[i]] <- WT_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>%
    as.data.frame()
  
  MUT_peaks[[i]] <- MUT_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>% 
    as.data.frame()
  
}

WT_gene_count <- list()

for (i in names(WT_peaks)) {
  WT_gene_count[[i]] <- WT_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}

MUT_gene_count <- list()

for (i in names(WT_peaks)) {
  MUT_gene_count[[i]] <- MUT_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}

common_gene_count <- list()

for (i in names(WT_peaks)) {
  common_gene_count[[i]] <- common_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}


all_counts <- list()


for (i in names(WT_gene_count)) {
  all_counts[[i]] <- full_join(full_join(WT_gene_count[[i]], MUT_gene_count[[i]], by = "SYMBOL"), common_gene_count[[i]], by = "SYMBOL")
  all_counts[[i]][is.na(all_counts[[i]])] <- 0 
  names(all_counts[[i]])[2:4] <- c("WT_count", "MUT_count", "Common_count")
}


for (i in names(all_counts)) {
  all_counts[[i]]$Unique_WT <- all_counts[[i]]$WT_count - all_counts[[i]]$Common_count
  all_counts[[i]]$Unique_MUT <- all_counts[[i]]$MUT_count - all_counts[[i]]$Common_count
  all_counts[[i]] <- filter(all_counts[[i]], all_counts[[i]]$Unique_MUT >= 0)
}


for (i in 1:length(all_counts)) {
  all_counts[[i]]$Perc_of_common <-  (all_counts[[i]]$Common_count/(all_counts[[i]]$Unique_WT + all_counts[[i]]$Unique_MUT + all_counts[[i]]$Common_count))*100
  all_counts[[i]]$Condition <- names(all_counts)[[i]]
}

for (i in 1:length(all_counts)) {
  all_counts[[i]] <-  filter(all_counts[[i]], all_counts[[i]]$Perc_of_common < 100 | all_counts[[i]]$Perc_of_common > 0 )
}




all_plots <- list() #plot data for all clusters

for (i in 1:10) {
  p <- ggplot(all_counts[[i]], aes(x=all_counts[[i]]$Unique_WT, y=all_counts[[i]]$Unique_MUT, color = all_counts[[i]]$Perc_of_common)) +
    geom_point() +
    scale_color_viridis()+
    geom_jitter(position = "jitter", aes(x=all_counts[[i]]$Unique_WT + 0.2, y=all_counts[[i]]$Unique_MUT + 0.2))+
    geom_smooth(method = "loess")+
    theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
          axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 15),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
          axis.title.y = element_text(face = "bold", color = "black", size = 15, vjust = 0),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          legend.position="top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Unique WT", y = "Unique MUT")+
    ggtitle(paste(names(all_counts)[i]))+
    xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))
  ggsave(p, filename = paste(names(all_counts)[i], ".png", sep = ""), width = 10, height = 10)
  
}



to_vio <- bind_rows(all_counts)
ggplot(filter(to_vio, to_vio$Common_count > 0), aes(x=Condition, y=Common_count)) + 
  geom_violin()+
  #geom_density()
  ylim(c(0, 10))

for (i in names(all_counts)) {
  density(all_counts[[i]]$MUT_count)
}

library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

DEGS <- read_excel_allsheets("DEGS_intra_cluster.xlsx")


for (i in intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- dplyr::left_join(all_counts[[i]], DEGS[[i]], by="SYMBOL")
}

for (i in intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- all_counts[[i]] %>%
    mutate(DEGS = case_when(
      avg_log2FC > 0 ~ "UP",
      avg_log2FC < 0 ~ "DOWN",
      is.na(avg_log2FC) ~ "NO_DEGS"))
  
}


for (i in intersect(names(DEGS), names(all_counts))) {
  p <- ggplot(all_counts[[i]], aes(x=all_counts[[i]]$Unique_WT, y=all_counts[[i]]$Unique_MUT, color = as.factor(DEGS))) +
    theme_classic()+
    geom_point() +
    geom_jitter(position = "jitter", aes(x=all_counts[[i]]$Unique_WT + 0.2, y=all_counts[[i]]$Unique_MUT + 0.2)) +
    labs(x = "Unique WT", y = "Unique MUT")+
    ggtitle(paste(i))+
    xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    scale_colour_manual(name = "Special Points",
                        values = c("blue", alpha("grey", 0.1), "red"))
  ggsave(p, filename = paste(i, ".png", sep = ""), width = 10, height = 10)
}



ggplot(all_counts[[1]], aes(x=all_counts[[1]]$Unique_WT, y=all_counts[[1]]$Unique_MUT, color = as.factor(DEGS))) +
  theme_classic()+
  geom_point() +
  labs(x = "Unique WT", y = "Unique MUT")+
  geom_jitter(position = "jitter", aes(x=all_counts[[1]]$Unique_WT + 0.2, y=all_counts[[1]]$Unique_MUT + 0.2))+
  xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
  ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
  scale_colour_manual(name = "Special Points",
                      values = c("blue", alpha("grey", 0.1), "red"))


for (i in  intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- all_counts[[i]][order(all_counts[[i]]$DEGS),]
}

to_write <- list()

for (i in intersect(names(DEGS), names(all_counts))) {
  df <- all_counts[[i]] %>%
    filter(DEGS!="NO_DEGS")
  to_write[[i]] <- df
}

openxlsx::write.xlsx(to_write, "Peaks_common_by_gene.xlsx", overwrite = F)


#Fig.6 k pseudotime during neural differentiation 