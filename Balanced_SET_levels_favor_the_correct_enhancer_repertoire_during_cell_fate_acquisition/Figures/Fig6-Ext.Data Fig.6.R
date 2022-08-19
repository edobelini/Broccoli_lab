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
Multi.object$Genotype <- str_split_fixed(Multi.object@meta.data$Condition, "_", 2)[,2]

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

mycol <- kelly(11)
mycol[1] <- "brown"

p1 = DimPlot(Multi.object, reduction = "umap", label = F, cols = mycol, pt.size = 1.5) +
  theme(plot.title = element_text(color = "black", size = 10, face = "bold.italic"),
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
p_val_adj_T = 0.05

cluster.markers <- FindAllMarkers(Multi.object,
                                 thresh.use = thresh.use,
                                 test.use = test.use,
                                 min.pct = min.pct,
                                 min.diff.pct = min.diff.pct,
                                 only.pos = LogFC.onlypos) %>%
                    as.data.frame() %>%
                    filter(p_val_adj < p_val_adj_T)

top100 <- cluster.markers %>%
group_by(cluster) %>%
top_n(n = 100, wt = avg_log2FC) %>%
as.data.frame()
rownames(top100) <- top100$gene

MG_list <- list()

for (i in 1:length(unique(Idents(Multi.object)))) {
  MG_list[[i]] <- filter(top100, top100$cluster == i)
  names(MG_list)[i] <- paste("Cluster", i, sep = "_")
}

openxlsx::write.xlsx(MG_list, "Marker_Genes_100.xlsx")

new.cluster.name <- c("ExN_UL", "OPCs_Oligo", "Astro", "IN", "Astro", "INP", "AP_RGC", "Late_Prog", "ExN_DL", 
  "Microglia", "ExN_UL", "Endothelium", "Endothelium", "Astro", "ExN_DL", "Late_Prog")



DEGS <- list()
Idents(Multi.object) <- Multi.object$Genotype

for (i in as.vector(unique(Multi.object$Label_cluster))) {
  filt <-  subset(Multi.object, subset = Label_cluster %in% i)
  DGE <- FindMarkers(filt, ident.1 = "mut" , ident.2 = "ctrl", 
                     logfc.threshold = 0.25, test.use = "wilcox")
  DEGS[[i]] <- as.data.frame(DGE)
  DEGS[[i]]$SYMBOL <- rownames(DEGS[[i]])
  names(DEGS)[i] <- paste("DEGS", i, sep = "_")
}

openxlsx::write.xlsx(MG_list, "DEGS_intra_cluster.xlsx")

A <- as.data.frame(prop.table(table(Idents(Multi.object), Multi.object$Condition), margin = 2))
A2 <- as.data.frame(table(Idents(Multi.object), Multi.object$Condition), margin = 2)
A$N <- A2$Freq


M <- ggplot(data = A, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity",  color = "black", position = "fill") +
  theme(plot.title = element_text(color = "black", size = 16, face = "bold.italic"),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size = 12, hjust = 1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "left",
        panel.background = element_rect(fill = "transparent", colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "fraction of cells")+
  scale_fill_manual(values = col)

png("plots/Perc_Condition_in_clusters_A.png",
    width     = 6,
    height    = 8,
    units     = "in",
    res       = 330)

M
dev.off()


library(monocle3)
########Trjectories with monocle3#######

cds_sort <- subset(Multi.object, subset=Condition %in% c("embryo_ctrl", "P2_ctrl") & Label_cluster %in% c("AP_RGC", "INP", "ExN_DL", "ExN_UL", "Late_Prog"))
cds <- as.cell_data_set(cds_sort, cell_metadata=as.data.frame(cds_sort@meta.data))
cds <- preprocess_cds(cds, num_dim = 10, norm_method = "log")
cds <- align_cds(cds, alignment_group = "Condition")

cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP")


cds@clusters$UMAP$partition <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$partition) <- rownames(cds_sort@meta.data)
cds@clusters$UMAP$clusters <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$clusters) <- rownames(cds_sort@meta.data)

traj_ctrl_Cond <- plot_cells(cds, color_cells_by = "Condition", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("ctrl")
traj_ctrl_Clus <- plot_cells(cds, color_cells_by = "cluster", group_cells_by = "cluster", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("ctrl")

Multi.object <- AddMetaData(object = Multi.object, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Pseudotime_ctrl" )

cds_sort <- subset(Multi.object, subset=Condition %in% c("embryo_mut", "P2_mut") & Label_cluster %in% c("AP_RGC", "INP", "ExN_DL", "ExN_UL", "Late_Prog"))
cds <- as.cell_data_set(cds_sort, cell_metadata=as.data.frame(cds_sort@meta.data))
cds <- preprocess_cds(cds, num_dim = 10, norm_method = "log")
cds <- align_cds(cds, alignment_group = "Condition")

cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP")


cds@clusters$UMAP$partition <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$partition) <- rownames(cds_sort@meta.data)
cds@clusters$UMAP$clusters <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$clusters) <- rownames(cds_sort@meta.data)


traj_ctrl_Cond <- plot_cells(cds, color_cells_by = "Condition", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("mut")
traj_ctrl_Clus <- plot_cells(cds, color_cells_by = "cluster", group_cells_by = "cluster", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("mut")

Multi.object <- AddMetaData(object = Multi.object, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Pseudotime_mut" )



library(scales)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(easyGgplot2)

ctrl_Time <- Multi.object@metadata %>% filter(!is.na(Pseudotime_ctrl))
mut_Time <- Multi.object@metadata %>% filter(!is.na(Pseudotime_mut))

ctrl_Time$Pseudotime_ctrl_scaled <- rescale(ctrl_Time$Pseudotime_ctrl, to = c(0, 100))
mut_Time$Pseudotime_mut_scaled <- rescale(mut_Time$Pseudotime_mut, to = c(0, 100))
ctrl_Time <- ctrl_Time[c(2, 4)]
ctrl_Time$Condition <- "ctrl"
names(ctrl_Time)[2] <- "Pseudotime"
mut_Time <- mut_Time[c(2, 4)]
mut_Time$Condition <- "mut"
names(mut_Time)[2] <- "Pseudotime"

pseudo <- rbind(ctrl_Time, mut_Time)
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_UL"] <- "ExN"
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_DL"] <- "ExN"
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_L1"] <- "ExN"

a <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("ExN", "INP" ,"AP_RGC")), 
                     xName='Pseudotime', groupName='Label_cluster', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = c("#BE0032","cornflowerblue" ,"#F3C300"))+
  #scale_fill_manual(values=kelly(length(unique(pseudo$Label_cluster))+1)[2:4])+
  theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "% of Cell")



b <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("ExN", "INP" ,"AP_RGC")), 
                     xName='Pseudotime', groupName='Condition', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = okabe(8)[6:7])+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Pseudotime", y = "% of Cell")



c <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("INP", "ExN", "Late_Prog", "AP_RGC", "IN", "ExN_L1", "Astro")),
                     xName='Pseudotime', groupName='Label_cluster', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = watlington(length(unique(pseudo$Label_cluster))))+
  #scale_fill_manual(values=kelly(length(unique(pseudo$Label_cluster))+1)[2:4])+
  theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "% of Cell")+
  ggtitle("Pseudotime distribution of lineage ExN -> INP -> AP_RGC")


png("ExN_INP_AP_RGC_in_pseudotime.png", res = 330, height = 20, width = 30, units = "in")
c/a/b
dev.off()

########################Diffrential Enhancer Usage#############

library(plyranges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(stringr)
library(viridis)
library(ggplot2)
library(openxlsx)
library(pheatmap)



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




all_plots <- list()

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





