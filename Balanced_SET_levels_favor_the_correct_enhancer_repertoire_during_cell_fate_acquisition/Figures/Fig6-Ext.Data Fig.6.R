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

#Fig.6 RNA velocity UMAP and heatmap 


