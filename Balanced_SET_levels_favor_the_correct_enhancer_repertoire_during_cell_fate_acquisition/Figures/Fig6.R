library(Seurat)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2)
library(openxlsx)
library(pals)
library(RColorBrewer)
library(openxlsx)
library(gprofiler2)
library(data.table)
library(readxl)
library(ComplexHeatmap)

options(future.globals.maxSize = 8000 * 1024^2)

setwd("./")

data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")

Multi.object <- CreateSeuratObject(counts = data$`Gene Expression`, project = "RNA", min.cells = 3, min.features = 200)
rm(data)

Multi.object$Condition <- as.data.frame(str_split_fixed(rownames(Multi.object@meta.data), pattern = "-", n = 2))[,2]
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "1"] <- "embryo_mut"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "2"] <- "embryo_ctrl"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "2"] <- "P2_ctrl"
Multi.object@meta.data["Condition"][Multi.object@meta.data["Condition"] == "2"] <- "P2_mut"

dir.create("plots/")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Multi.object[["percent.mt"]] <- PercentageFeatureSet(Multi.object, pattern = "^mt-")
VlnPlot(Multi.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = F, group.by = "Condition")+
  ggsave("plots/vio_plot.png", dpi = 330, scale = 0.8, width = 12, height = 5)

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
p1

p2 = DimPlot(Multi.object , reduction = "umap", group.by = "Condition", label = F, pt.size = 0.1, cols = c("red", "blue", "green", "black"))+
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=20, hjust =1),
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        legend.position="left",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")
p2

png("plots/UMAP.png", width = 5, height = 10, units = "cm", res = 330, pointsize = 3)
p1+p2
dev.off()

bench_in <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
bench_plot_list <- list()
for (i in bench_in) {
  FeaturePlot(Multi.object, features = c(i))+
    scale_colour_gradientn(colours = viridis(20))
  ggsave(paste("plots/",i ,".png"), dpi = 330, scale = 0.8, width = 12, height = 5)
}

thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
test.use = "wilcox"
LogFC.onlypos = TRUE
# find markers
cluster.markers = FindAllMarkers(Multi.object,
                                 thresh.use = thresh.use,
                                 test.use=test.use,
                                 min.pct=min.pct,
                                 min.diff.pct=min.diff.pct,
                                 only.pos=LogFC.onlypos)

top100 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()
rownames(top100) <- top100$gene
write.csv(cluster.markers, "Marker_genes.csv", quote = FALSE)

MG_list <- list()

for (i in 1:length(unique(Idents(Multi.object)))) {
  MG_list[[i]] <- filter(top100, top100$cluster==i)
  names(MG_list)[i] <- paste("Cluster", i, sep = "_")
}

openxlsx::write.xlsx(MG_list, "Marker_Genes_100.xlsx",asTable = F)
