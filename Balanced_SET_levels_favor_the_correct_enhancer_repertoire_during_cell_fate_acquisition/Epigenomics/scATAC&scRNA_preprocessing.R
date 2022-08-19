library(ArchR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plyranges)
library(Seurat)

#set Number of threads
addArchRThreads(threads = 30)

#Reference genome

addArchRGenome("mm10")

#Load scRNA object and Seurat object creation 

data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")

Multi.object <- CreateSeuratObject(counts = data$`Gene Expression`, project = "RNA", min.cells = 3, min.features = 200)
rm(data)

Multi.object$Condition <- as.data.frame(str_split_fixed(rownames(Multi.object@meta.data), pattern = "-", n = 2))[,2]
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "1"] <- "embryo_mut"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "2"] <- "embryo_ctrl"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "3"] <- "P2_ctrl"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "4"] <- "P2_mut"

dir.create("plots/")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Multi.object[["percent.mt"]] <- PercentageFeatureSet(Multi.object, pattern = "^mt-")
VlnPlot(Multi.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = F, group.by = "Condition")+
  ggsave("plots/vio_plot.png", dpi = 330, scale = 0.8, width = 12, height = 5)

min_nFeature_RNA = 200
max_nFeature_RNA = 4000
max_percent_MT = 10

  # calculate %mt features
Multi.object[["percent.mt"]] <- PercentageFeatureSet(Multi.object, pattern = "^mt-")
  # show violin plots to set quality threshold
print(VlnPlot(Multi.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 1))
  # subset cells matching quality requirment
Multi.object <-  subset(Multi.object, subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)
  # Violin plots after filtering
print(VlnPlot(Multi.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 1))
  # Data normalization
Multi.object <- NormalizeData(Multi.object, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Multi.object)
Multi.object <- ScaleData(Multi.object, features = all.genes)


Multi.object <- RunPCA(Multi.object, features = VariableFeatures(object = Multi.object)) %>% 
    JackStraw(Multi.object, num.replicate = 100) %>% 
    ScoreJackStraw(Multi.object, dims = 1:10) %>%
    ElbowPlot(Multi.object) %>% 
    FindNeighbors(Multi.object, dims = 1:20) %>%
    FindClusters(Multi.object, resolution = 0.5) %>%
    RunUMAP(Multi.object, dims = 1:15)

colors <- kelly(11)
colors[1] <- "brown" 

p1 = DimPlot(Multi.object, reduction = "umap", label = F, cols = mycol, pt.size = 1.5) +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=20, hjust =1),
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("UMAP 1:15, FindNeighbors 1:20")




#scATAC data preprocessing 

#create Arrow files for ArchR preprocessing 

ArrowFiles <- createArrowFiles(
  inputFiles = c("Documents/scATAC_Multiome_Bam/Embryo_Ctrl/atac_fragments.tsv.gz",
                 "Documents/scATAC_Multiome_Bam/Embryo_Mut/atac_fragments.tsv.gz",
                 "Documents/scATAC_Multiome_Bam/P2_Ctrl/atac_fragments.tsv.gz",
                 "Documents/scATAC_Multiome_Bam/P2_Mut/atac_fragments.tsv.gz"),
  sampleNames = c("embryo_ctrl","embryo_mut","P2_ctrl","P2_mut"),
  minTSS = 0, #We do not use filter here becasue we filter based on scRNA cells
  minFrags = 1, #We do not use filter here becasue we filter based on scRNA cells 
  maxFrags = 10000000, #We do not use filter here becasue we filter based on scRNA cells
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


#Create ArrowProject



scATAC_Multiome <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Documents/scATAC_Multiome_Bam/Table_heatmap/",
  copyArrows = F 
)

#add RNA matrix to ArchRProject

seRNA <- import10xFeatureMatrix(
  input = c("Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Multiome_SETBP1/ATAC_CRE_Plus_embrioni/Embryo_Mut/outs/filtered_feature_bc_matrix.h5",
            "Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Multiome_SETBP1/ATAC__CRE_Minus_embrioni/Embryo_Ctrl/outs/filtered_feature_bc_matrix.h5",
            "Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Multiome_SETBP1/ATAC_ctrl_neonati/P2_Ctrl/outs/filtered_feature_bc_matrix.h5",
            "Share_HSR/Ric.Broccoli/zaghi.mattia/SETBP1_epigenomics/Multiome_SETBP1/ATAC_mutato_neonati/P2_Mut/outs/filtered_feature_bc_matrix.h5"),
  names = c("embryo_ctrl","embryo_mut","P2_ctrl","P2_mut"))




seRNAcombined<-cbind(assay(seRNA[[1]]), assay(seRNA[[2]]), assay(seRNA[[3]]), assay(seRNA[[4]]))

seRNA <- SummarizedExperiment(assays=list(counts=seRNAcombined), rowRanges= rowRanges(seRNA[[1]])) 
