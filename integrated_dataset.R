library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(scCustomize)
set.seed(1234)
gc()

data1.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5833_2_output_filtered.h5')
data2.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5807_3_output_filtered.h5') 
data3.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5285_2_output_filtered.h5') 
data4.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5880_1_output_filtered.h5') 
data5.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5880_2_output_filtered.h5') 
data6.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5833_1_output_filtered.h5')
data7.data <- Read10X_h5('~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/5285_1_output_filtered.h5') 


data1 <- CreateSeuratObject(counts = data1.data, project = 'Unchallenged')
data2 <- CreateSeuratObject(counts = data2.data, project = 'SD Bleo')
data3 <- CreateSeuratObject(counts = data3.data, project = 'Rep Bleo')
data4 <- CreateSeuratObject(counts = data4.data, project = 'Sftpc-trace')
data5 <- CreateSeuratObject(counts = data5.data, project = 'Scgb1a1-trace')
data6 <- CreateSeuratObject(counts = data6.data, project = 'Hif1/2 Unchallenged')
data7 <- CreateSeuratObject(counts = data7.data, project = 'Hif1/2 Rep Bleo')

adata <- merge(data1, y=c(data2,data3,data4, data5,data6,data7))

adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^mt-")
adata[["percent.ribo"]] <- PercentageFeatureSet(adata, pattern = "^Rp")
VlnPlot(adata, features = c('percent.mt', 'percent.ribo'))

adata <- subset(adata, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
adata <- SCTransform(adata, vars.to.regress = c('percent.mt', 'percent.ribo'), return.only.var.genes = F)
#adata <- NormalizeData(adata)
#adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
#adata <- ScaleData(adata, vars.to.regress = c('percent.mt', 'percent.ribo', 'nCount_RNA'), n_jobs=8)
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DotPlot_scCustom(adata, features = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Ms4a1', 'Mki67', 'Krt5', 'Sox8'))

#remove non-epithelial
adata <- subset(adata, idents = c(0:14,17:22))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')

DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')

epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Cdc20', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
adata <- subset(adata, idents = c(0:4, 6:11,13,15:20))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')

DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')

epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Cdc20', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets/low quality
adata <- subset(adata, idents = c(0:17))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')

DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')

epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Cdc20', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets/low quality
adata <- subset(adata, idents = c(0:11,13:16))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')

DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')

epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Cdc20', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))

DefaultAssay(adata) <- 'RNA'
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata, vars.to.regress = c('percent.mt', 'percent.ribo', 'nCount_RNA'), n_jobs=8)
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')
epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', 'Cdc20', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lyz1', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Vim', 'Fn1', 'Itgav', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

#remove low quality 6,8, doublet 18
adata <- subset(adata, idents = c(0:5,7,9:17))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')
epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', 'Cdc20', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lyz1', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Vim', 'Fn1', 'Itgav', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

#remove low quality 5
adata <- subset(adata, idents = c(0:4,6:16))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')
epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', 'Cdc20', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lyz1', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Vim', 'Fn1', 'Itgav', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

#remove low quality 13
adata <- subset(adata, idents = c(0:12,14:15,17,18))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')
epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', 'Cdc20', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lyz1', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Vim', 'Fn1', 'Itgav', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

#remove low quality 11
adata <- subset(adata, idents = c(0:11,13,14))
adata <- RunPCA(adata, verbose = F)
adata <- RunHarmony(adata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata <- FindNeighbors(adata, reduction = 'harmony', dims = 1:45)
adata <- FindClusters(adata, resolution = 0.9)
adata <- RunUMAP(adata, dims = 1:45, reduction = 'harmony')


DimPlot_scCustom(adata, figure_plot = T)
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')
epi_genes = c('Ptprc', 'Col1a1', 'Pecam1', 'Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', 'Cdc20', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lyz1', 'H2-Aa', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Gdf15', 'Hopx', 'Rtkn2', 'Ager', 'Vegfa', 'Col4a3', 'Calca', 'Mki67', 'Vim', 'Fn1', 'Itgav', 'Sox8')
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

Idents(adata) <- 'seurat_clusters'
adata <- FindSubCluster(adata, 4, graph.name = 'RNA_snn', resolution = 0.3, subcluster.name='subcluster')
Idents(adata) <- 'subcluster'
adata <- FindSubCluster(adata, 10, graph.name = 'RNA_snn', resolution = 0.3, subcluster.name='subcluster1')
Idents(adata) <- 'subcluster1'
adata <- FindSubCluster(adata, 7, graph.name = 'RNA_snn', resolution = 0.3, subcluster.name='subcluster2')
Idents(adata) <- 'subcluster2'
adata <- FindSubCluster(adata, 11, graph.name = 'RNA_snn', resolution = 0.3, subcluster.name='subcluster3')
Idents(adata) <- 'subcluster3'

DimPlot_scCustom(adata, figure_plot = T)
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')


Idents(adata) <- 'subcluster3'
Idents(adata, cells = WhichCells(adata, idents = c('10_3'))) <- "M-like cells"
Idents(adata, cells = WhichCells(adata, idents = c('10_1'))) <- "Secretory"
Idents(adata, cells = WhichCells(adata, idents = c('10_0'))) <- "Mucous secretory"
Idents(adata, cells = WhichCells(adata, idents = c(12))) <- "Proliferating"
Idents(adata, cells = WhichCells(adata, idents = c(13))) <- "PNEC"
Idents(adata, cells = WhichCells(adata, idents = c(2,5))) <- "MCC"
Idents(adata, cells = WhichCells(adata, idents = c('10_2'))) <- "Airway intermediate"
Idents(adata, cells = WhichCells(adata, idents = c('4_0'))) <- "Alveolar intermediate"
Idents(adata, cells = WhichCells(adata, idents = c('4_2'))) <- "Early intermediate"
Idents(adata, cells = WhichCells(adata, idents = c('11_0', '11_1'))) <- "Basal-like"
Idents(adata, cells = WhichCells(adata, idents = c('7_0', '7_1', '7_2'))) <- "BASC-like"
Idents(adata, cells = WhichCells(adata, idents = c(0,1,3,6,8,9))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c('4_1'))) <- "AT1"
adata$celltype <- Idents(adata)

Idents(adata) <- 'subcluster3'
Idents(adata, cells = WhichCells(adata, idents = c('4_1'))) <- "AT1"
Idents(adata, cells = WhichCells(adata, idents = c(0,1,3,6,8,9))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c('11_0', '11_1'))) <- "Basal-like"
Idents(adata, cells = WhichCells(adata, idents = c('4_2'))) <- "Early intermediate"
Idents(adata, cells = WhichCells(adata, idents = c('4_0'))) <- "Alveolar intermediate"
Idents(adata, cells = WhichCells(adata, idents = c('10_2'))) <- "Airway intermediate"
Idents(adata, cells = WhichCells(adata, idents = c(2,5))) <- "MCC"
Idents(adata, cells = WhichCells(adata, idents = c(13))) <- "PNEC"
Idents(adata, cells = WhichCells(adata, idents = c(12))) <- "Proliferating"
Idents(adata, cells = WhichCells(adata, idents = c('10_0'))) <- "Mucous secretory"
Idents(adata, cells = WhichCells(adata, idents = c('10_1'))) <- "Secretory"
Idents(adata, cells = WhichCells(adata, idents = c('10_3'))) <- "M-like cells"
adata$celltype_rev <- Idents(adata)

Idents(adata) <- 'celltype'

DimPlot_scCustom(adata, figure_plot = T)
DotPlot_scCustom(adata, features = epi_genes, cluster.idents=TRUE)+theme(axis.text.x = element_text(angle = 45, hjust=1))
VlnPlot(adata, features = 'nFeature_RNA')

Idents(adata) <- 'celltype'
saveRDS(adata, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/total_dataset_20230718.rds')


Idents(adata) <- 'orig.ident'
fig1 <- subset(adata, idents = c('Unchallenged', 'SD Bleo', 'Rep Bleo'))
fig1 <- RunPCA(fig1, verbose = F)
fig1 <- RunHarmony(fig1, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig1 <- FindNeighbors(fig1, reduction = 'harmony', dims = 1:45)
fig1 <- RunUMAP(fig1, dims = 1:45, reduction = 'harmony')
Idents(fig1) <- 'celltype'
DimPlot_scCustom(fig1, figure_plot = T)

fig2 <- subset(adata, idents = c('Unchallenged', 'Hif1/2 Unchallenged', 'Rep Bleo', 'Hif1/2 Rep Bleo'))
fig2 <- RunPCA(fig2, verbose = F)
fig2 <- RunHarmony(fig2, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig2 <- FindNeighbors(fig2, reduction = 'harmony', dims = 1:45)
fig2 <- RunUMAP(fig2, dims = 1:45, reduction = 'harmony')
Idents(fig2) <- 'celltype'
DimPlot_scCustom(fig2, figure_plot = T)


fig3 <- subset(adata, idents = c('Rep Bleo', 'Sftpc-trace', 'Scgb1a1-trace'))
fig3 <- RunPCA(fig3, verbose = F)
fig3 <- RunHarmony(fig3, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig3 <- FindNeighbors(fig3, reduction = 'harmony', dims = 1:45)
fig3 <- RunUMAP(fig3, dims = 1:45, reduction = 'harmony')
Idents(fig3) <- 'celltype'
DimPlot_scCustom(fig3, figure_plot = T)


saveRDS(fig1, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/wtsdrep_20230718.rds')
saveRDS(fig2, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/hif_20230718.rds')
saveRDS(fig3, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/trace_20230718.rds')

######################################################have not run yet#######
Idents(fig1) <- 'celltype'

ct_prop_table <- as.data.frame(prop.table(table(Idents(fig1), fig1$orig.ident), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )
ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Unchallenged', 'SD Bleo', 'Rep Bleo'))



p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")

p1 + scale_fill_manual(values = glasbey_pal)


#fig2
DimPlot_scCustom(adata, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 28, palette = "glasbey"), pt.size = 0.1) + theme(plot.title = element_blank())
DimPlot_scCustom(adata, figure_plot = T, group.by = 'orig.ident')


ct_prop_table <- as.data.frame(prop.table(table(Idents(fig2), fig2$orig.ident), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )

ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Unchallenged', 'Hif1/2 Unchallenged', 'Rep Bleo', 'Hif1/2 Rep Bleo'))


p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")
p1 + scale_fill_manual(values = glasbey_pal)

#fig2
cdata <- RunPCA(cdata, verbose = F)
cdata <- RunHarmony(cdata, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
cdata <- FindNeighbors(cdata, reduction = 'harmony', dims = 1:45)
cdata <- RunUMAP(cdata, dims = 1:45, reduction = 'harmony')
Idents(cdata) <- 'celltype'
DimPlot_scCustom(cdata, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 28, palette = "glasbey"), pt.size = 0.1) + theme(plot.title = element_blank())
saveRDS(cdata, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/fig2_20230710.rds' )
ct_prop_table <- as.data.frame(prop.table(table(Idents(cdata), cdata$orig.ident), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )
ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Unchallenged', 'Hif1/2 Unchallenged', 'Rep Bleo', 'Hif1/2 Rep Bleo'))


p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

p1 + scale_fill_manual(values = glasbey_pal)

#fig3
trace <- RunPCA(trace, verbose = F)
trace <- RunHarmony(trace, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
trace <- FindNeighbors(trace, reduction = 'harmony', dims = 1:45)
trace <- RunUMAP(trace, dims = 1:45, reduction = 'harmony')
Idents(trace) <- 'celltype'
DimPlot_scCustom(trace, label = F, figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 28, palette = "glasbey"), pt.size = 0.1) + theme(plot.title = element_blank())
saveRDS(trace, file = '~/Kropski_Lab Dropbox/Jon Kropski/Lab Data/hif_paper_10x/cb_outs/trace_20230710.rds' )
ct_prop_table <- as.data.frame(prop.table(table(Idents(trace), trace$orig.ident), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )
ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Sftpc-trace', 'Scgb1a1-trace'))


p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")
p1 + scale_fill_manual(values = glasbey_pal)
