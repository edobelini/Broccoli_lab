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


