library(Seurat)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2)
library(openxlsx)
library(pals)
library(RColorBrewer)


options(future.globals.maxSize = 8000 * 1024^2)

setwd("/home/edoardo/sno116/Multi/")
load("Multiome.RData")

#load("Multiome.RData")

pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

pbmc_RNA <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "RNA", min.cells = 3, min.features = 200)
#pbmc_ATAC <- CreateSeuratObject(counts = pbmc.data$Peaks, project = "ATAC", min.cells = 3, min.features = 200)
rm(pbmc.data)
pbmc_RNA$Condition <- as.data.frame(str_split_fixed(rownames(pbmc_RNA@meta.data), pattern = "-", n = 2))[,2]

pbmc_RNA@meta.data["Condition"][pbmc_RNA@meta.data["Condition"] == "1"] <- "Ube3a_mut"
pbmc_RNA@meta.data["Condition"][pbmc_RNA@meta.data["Condition"] == "2"] <- "WT"


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc_RNA[["percent.mt"]] <- PercentageFeatureSet(pbmc_RNA, pattern = "^mt-")
VlnPlot(pbmc_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = F, group.by = "Condition")+
  ggsave("plots/vio_plot.png", dpi = 330, scale = 0.8, width = 12, height = 5)

#pbmc_RNA <- subset(pbmc_RNA, subset = nFeature_RNA > 1000)

pbmc_RNA <- NormalizeData(pbmc_RNA, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc_RNA)
pbmc_RNA <- ScaleData(pbmc_RNA, features = all.genes)


pbmc_RNA <- RunPCA(pbmc_RNA, features = VariableFeatures(object = pbmc_RNA))

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc_RNA <- JackStraw(pbmc_RNA, num.replicate = 100)
pbmc_RNA <- ScoreJackStraw(pbmc_RNA, dims = 1:10)
ElbowPlot(pbmc_RNA)

pbmc_RNA <- FindNeighbors(pbmc_RNA, dims = 1:20)
pbmc_RNA <- FindClusters(pbmc_RNA, resolution = 0.5)

mycol <- as.character(polychrome(length(unique(Idents(pbmc_RNA)))))

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc_RNA <- RunUMAP(pbmc_RNA, dims = 1:15)
View(pbmc_RNA@meta.data)
Idents(pbmc_RNA) <- pbmc_RNA$seurat_clusters
pbmc_RNA@meta.data$nCount_RNA <- pbmc_RNA$nCount_RNA
pbmc_RNA@meta.data$nFeature_RNA <- pbmc_RNA$nFeature_RNA

p1 = DimPlot(pbmc_RNA, reduction = "umap", label = T, pt.size = 0.1, cols = mycol) +
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

p2 = DimPlot(pbmc_RNA , reduction = "umap", group.by = "Condition",label = F, pt.size = 0.1, cols = c("red", "blue", "green", "black"))+
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

png("plots/UMAP.png", width = 12, res = 330, height = 5)
p1+p2
dev.off()

bench_in <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
bench_plot_list <- list()
for (i in bench_in) {
  FeaturePlot(pbmc_RNA, features = c(i))+
    scale_colour_gradientn(colours = viridis(20))
  ggsave(paste("plots/",i ,".png"), dpi = 330, scale = 0.8, width = 12, height = 5)
}

thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
test.use = "wilcox"
LogFC.onlypos = TRUE
# find markers
cluster.markers = FindAllMarkers(pbmc_RNA,
                                 thresh.use = thresh.use,
                                 test.use=test.use,
                                 min.pct=min.pct,
                                 min.diff.pct=min.diff.pct,
                                 only.pos=LogFC.onlypos)

top100 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()
rownames(top100) <- top100$gene
write.csv(cluster.markers, "Marker_genes.csv", quote = FALSE)

MG_list <- list()

for (i in 1:length(unique(Idents(pbmc_RNA)))) {
  MG_list[[i]] <- filter(top100, top100$cluster==i)
  names(MG_list)[i] <- paste("Cluster", i, sep = "_")
}

openxlsx::write.xlsx(MG_list, "Marker_Genes_100.xlsx",asTable = F)


