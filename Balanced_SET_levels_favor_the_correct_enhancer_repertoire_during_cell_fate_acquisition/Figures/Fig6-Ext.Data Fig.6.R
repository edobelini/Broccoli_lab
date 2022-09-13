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


#Fig.6 h Expression plots during pseudotime of genes 

