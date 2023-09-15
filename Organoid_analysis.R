library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(Signac)
library(harmony)
library(tidyverse)
library(patchwork)
library(viridis)
library(scCustomize)
library(qs)
library(dplyr)
library(MuDataSeurat)
set.seed(1234)

data1.data <- Read10X_h5('~/Desktop/organoid/filtered_feature_bc_matrix.h5')
data1 <- CreateSeuratObject(counts = data1.data, project = "Vehicle", min.cells = 3, min.features = 500)

data2.data <- Read10X_h5('~/Desktop/organoid/filtered_feature_bc_matrix2.h5')
data2 <- CreateSeuratObject(counts = data2.data, project = "PT-2385", min.cells = 3, min.features = 500)

adata <- merge(data1, y=data2)

adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
adata[["percent.ribo"]] <- PercentageFeatureSet(adata, pattern = "^RP")
VlnPlot(adata, features = c('percent.mt', 'percent.ribo'))

DefaultAssay(adata) <- 'RNA'
adata <- subset(adata, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata, vars.to.regress = c('percent.mt', 'percent.ribo', 'nCount_RNA'), n_jobs=8)
adata <- RunPCA(adata, verbose = F)
#adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 12)
adata <- FindNeighbors(adata,dims = 1:18)
adata <- FindClusters(adata, resolution = 0.6)
adata <- RunUMAP(adata, dims = 1:18)
DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')

VlnPlot_scCustom(adata, features = 'nFeature_RNA')

#remove cluster 7 low counts
adata <- subset(adata, idents = c(0:6,8))
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata, vars.to.regress = c('percent.mt', 'percent.ribo', 'nCount_RNA'), n_jobs=8)
adata <- RunPCA(adata, verbose = F)
#adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 12)
adata <- FindNeighbors(adata,dims = 1:12)
adata <- FindClusters(adata, resolution = 0.6)
adata <- RunUMAP(adata, dims = 1:12)
DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')


epi_genes4<-c('NKX2-1','CD24', 'TP63', 'RNASE1', 'MUC5B', 'SCGB3A1', 'SCGB3A2', 'CD74', 'HLA-DRA', 'PGC', 'GPC5', 'ETV5', 'LRRK2', 'NEDD4L', 'APOE', 'SLPI','SFTPA2', 'SFTPB', 'SFTPC', 'SFTPD', 'LAMP3', 'ABCA3', 'NAPSA','SLC34A2', 'FOLR1', 'CD63', 'MYL6', 'TNNC1', 'SLC7A11', 'ADARB2', 'ERBB4', 'CEACAM6', 'KRT8', 'CLDN3', 'CLDN4', 'LGALS3', 'LGALS1', 'CDKN1A', 'GDF15', 'AGER', 'HOPX','MKI67', 'FGFR2', 'FN1', 'COL1A1', 'TM4SF1', 'HIF1A', 'EPAS1','SLC2A1', 'CA9', 'VEGFA','SOX4','VIM', 'IGFBP2', 'IGFBP5')

DotPlot_scCustom(adata, features = epi_genes4, cluster.idents = T)+theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(adata) <- 'seurat_clusters'
Idents(adata, cells = WhichCells(adata, idents = c(7))) <- "Basal"
Idents(adata, cells = WhichCells(adata, idents = c(4))) <- "Aberrant intermediate"
Idents(adata, cells = WhichCells(adata, idents = c(5))) <- "AT2 intermediate"
Idents(adata, cells = WhichCells(adata, idents = c(6))) <- "Proliferating"
Idents(adata, cells = WhichCells(adata, idents = c(2,3))) <- "RASC-AT2"
Idents(adata, cells = WhichCells(adata, idents = c(1))) <- "AT2 - mature"
Idents(adata, cells = WhichCells(adata, idents = c(0))) <- "AT2 - activated"
adata$celltype <- Idents(adata)


Idents(adata) <- 'seurat_clusters'
Idents(adata, cells = WhichCells(adata, idents = c(0))) <- "AT2 - activated"
Idents(adata, cells = WhichCells(adata, idents = c(1))) <- "AT2 - mature"
Idents(adata, cells = WhichCells(adata, idents = c(2,3))) <- "RASC-AT2"
Idents(adata, cells = WhichCells(adata, idents = c(6))) <- "Proliferating"
Idents(adata, cells = WhichCells(adata, idents = c(5))) <- "AT2 intermediate"
Idents(adata, cells = WhichCells(adata, idents = c(4))) <- "Aberrant intermediate"
Idents(adata, cells = WhichCells(adata, idents = c(7))) <- "Basal"
adata$celltype_rev <- Idents(adata)

DotPlot_scCustom(adata, features = epi_genes4, cluster.idents = T)+theme(axis.text.x = element_text(angle = 45, hjust=1))
DimPlot_scCustom(adata, figure_plot = T)

epi_genes5 <- c('CD24', 'TP63', 'CLDN10', 'SCGB3A2', 'SFTPB', 'ETV5', 'LRRK2', 'FABP5', 'APOE', 'SFTPC', 'LAMP3', 'EPHA4', 'GPC5', 'MKI67',  'CLDN4', 'VIM', 'IGFBP2', 'FN1', 'LGALS3', 'CP', 'CXCL8')
DotPlot_scCustom(adata, features = epi_genes5, cluster.idents = F)+theme(axis.text.x = element_text(angle = 45, hjust=1))

markers <- FindAllMarkers(adata)
write.csv(markers, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/McCall/PT_organoid/analysis/pt_organoid_20230501_markers.csv')

saveRDS(adata, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/McCall/PT_organoid/analysis/pt_organoid_20230501.rds')
