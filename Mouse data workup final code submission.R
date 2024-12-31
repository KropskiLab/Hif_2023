library(Seurat)
library(harmony)
library(scCustomize)
set.seed(1234)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 40000 * 1024^2)
library(RColorBrewer)
library(viridis)
library(Nebulosa)
library(alluvial)
library(enrichR)
library(org.Mm.eg.db)

setwd("C:/Users/mccalas1/Kropski_Lab Dropbox/Scott McCall/cb_outs/20230718")

adata <- readRDS(file ="C:/Users/mccalas1/Kropski_Lab Dropbox/Scott McCall/cb_outs/20230718/total_dataset_20230718.rds")

#collapsing early and alveolar intermediates for simplicity and reordering, celltype2
Idents(adata) <- 'celltype'
Idents(adata, cells = WhichCells(adata, idents = c("M-like cells"))) <- "M-like cells"
Idents(adata, cells = WhichCells(adata, idents = c("Secretory"))) <- "Secretory"
Idents(adata, cells = WhichCells(adata, idents = c("Mucous secretory"))) <- "Secretory-Muc5b+"
Idents(adata, cells = WhichCells(adata, idents = c("BASC-like"))) <- "BASC-like"
Idents(adata, cells = WhichCells(adata, idents = c("Proliferating"))) <- "Proliferating"
Idents(adata, cells = WhichCells(adata, idents = c("PNEC"))) <- "PNEC"
Idents(adata, cells = WhichCells(adata, idents = c("MCC"))) <- "MCC"
Idents(adata, cells = WhichCells(adata, idents = c("Basal-like"))) <- "Basal-like"
Idents(adata, cells = WhichCells(adata, idents = c("Airway intermediate"))) <- "Airway intermediate"
Idents(adata, cells = WhichCells(adata, idents = c("Early intermediate", "Alveolar intermediate"))) <- "Alveolar intermediate"
Idents(adata, cells = WhichCells(adata, idents = c("AT2"))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c("AT1"))) <- "AT1"
adata$celltype2 <- Idents(adata)
Idents(adata) <- 'celltype2'
#adata$celltype <- factor(adata$celltype,levels=c("AT1","AT2","Alveolar intermediate", "Secretory", "Secretory-Muc5b+", "BASC-like", "Airway intermediate", "Basal-like", "MCC", "PNEC", "M-like cells", "Proliferating"))

Idents(adata) <- 'celltype'
DimPlot_scCustom(adata, figure_plot = T)

Idents(adata) <- 'celltype2'
DimPlot_scCustom(adata, figure_plot = T, pt.size = 3, label = F, raster = T)

Idents(adata) <- 'orig.ident'
DimPlot_scCustom(adata,  label = F, split.by = 'orig.ident', raster = TRUE, pt.size = 6.0, colors_use = c("#757576", "#C49A6C", "#010101","#92298D","#369846", "#cc7722","#0000FF"))
DimPlot_scCustom(adata,  label = F, figure_plot = TRUE , raster = TRUE, pt.size = 6.0, colors_use = c("#757576", "#C49A6C", "#010101","#92298D","#369846", "#cc7722","#0000FF"))
Idents(adata) <- 'celltype2'
col_var <- DiscretePalette_scCustomize(num_colors = 12, palette = "varibow")
DimPlot_scCustom(adata,  label = F, figure_plot = TRUE , raster = TRUE, pt.size = 6.0, colors_use = col_var)
Plot_Density_Custom(seurat_object = adata, features = c("Hif1a", "Epas1"))




epi_genes_figure = c('Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cyp2f2', 'Cd74','Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Krt8', 'Cdkn1a','Lgals3', 'Hopx', 'Rtkn2', 'Ager', 'Sox8' , 'Calca', 'Mki67')

DotPlot(adata, features = epi_genes_figure, cluster.idents= TRUE, dot.scale = 6) +
  RotatedAxis() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

ct_prop_table <- as.data.frame(prop.table(table(Idents(adata), adata$orig.ident), margin = 2))
ct_prop_table <- ct_prop_table %>% 
  rename(
    Cell_type = Var1,
    Condition = Var2
  )
ct_prop_table$Condition <- factor(ct_prop_table$Condition, levels = c('Unchallenged', 'SD Bleo', 'Rep Bleo', 'Scgb1a1-trace', 'Sftpc-trace', 'Hif1/2 Unchallenged', 'Hif1/2 Rep Bleo'))
write.csv(ct_prop_table,file = "revamp_prop_table_20230718.csv")

table(Idents(adata),adata$orig.ident)



p1<- ggplot(ct_prop_table, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Cell Types by Condition") +
  xlab("")

glasbey_pal <- DiscretePalette_scCustomize(num_colors = 32, palette = "glasbey")

p1 + scale_fill_manual(values = glasbey_pal)


epi_genes = c('Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cd74', 'Lrrk2', 'Lrp2', 'Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Abca3', 'Nckap5', 'Krt8', 'Cdkn1a', 'Cldn4', 'Hopx', 'Rtkn2', 'Ager', 'Calca', 'Mki67', 'Hif1a','Epas1')
epi_genes_figure = c('Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cyp2f2', 'Cd74','Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Krt8', 'Cdkn1a', 'Hopx', 'Rtkn2', 'Ager', 'Sox8' , 'Calca', 'Mki67')
write.csv(epi_genes, file = "epi_genes.csv")


Idents(adata) <- 'celltype2'
total_group_markers <-FindAllMarkers(adata, test.use = 'negbinom', logfc.threshold = 0.5)
write.csv(total_group_markers, file = "final_markers_total_set.csv")

Idents(adata) <- 'orig.ident'
fig1 <- subset(adata, idents = c('Unchallenged', 'SD Bleo', 'Rep Bleo'))
fig1 <- RunPCA(fig1, verbose = F)
fig1 <- RunHarmony(fig1, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig1 <- FindNeighbors(fig1, reduction = 'harmony', dims = 1:45)
fig1 <- RunUMAP(fig1, dims = 1:45, reduction = 'harmony')
Idents(fig1) <- 'celltype2'
DimPlot_scCustom(fig1, figure_plot = T)

fig2 <- subset(adata, idents = c('Unchallenged', 'Hif1/2 Unchallenged', 'Rep Bleo', 'Hif1/2 Rep Bleo'))
fig2 <- RunPCA(fig2, verbose = F)
fig2 <- RunHarmony(fig2, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig2 <- FindNeighbors(fig2, reduction = 'harmony', dims = 1:45)
fig2 <- RunUMAP(fig2, dims = 1:45, reduction = 'harmony')
Idents(fig2) <- 'celltype2'
DimPlot_scCustom(fig2, figure_plot = T)


fig3 <- subset(adata, idents = c('Rep Bleo', 'Sftpc-trace', 'Scgb1a1-trace'))
fig3 <- RunPCA(fig3, verbose = F)
fig3 <- RunHarmony(fig3, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig3 <- FindNeighbors(fig3, reduction = 'harmony', dims = 1:45)
fig3 <- RunUMAP(fig3, dims = 1:45, reduction = 'harmony')
Idents(fig3) <- 'celltype2'
DimPlot_scCustom(fig3, figure_plot = T)

saveRDS(fig1, file = "Fig1_20230727.rds")
saveRDS(fig2, file = "Fig2_20230727.rds")
saveRDS(fig3, file = "Fig3_20230727.rds")
saveRDS(adata, file = "totalMouse_20230727.rds")

#Figure 1 data workup for Dimplots, dotplots, and alluvial

Idents(fig1) <- 'orig.ident'
DimPlot_scCustom(fig1,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6.0, colors_use = c("#757576", "#C49A6C", "#010101"))

table(Idents(fig1),fig1$celltype2)


col_var <- DiscretePalette_scCustomize(num_colors = 12, palette = "varibow")
dat_col_var <- as.data.frame(col_var)
write.csv(dat_col_var, file = "varibowColorList.csv")

Idents(fig1) <- 'celltype2'
DimPlot_scCustom(fig1,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6.0, colors_use = col_var)

Idents(fig1) <- 'celltype2'
epi_genes_figure = c('Epcam', 'Nkx2-1', 'Sox2', 'Krt5', 'Trp63', 'Foxj1', "Scgb1a1", 'Scgb3a2', 'Muc5b', 'Cyp2f2', 'Cd74','Sftpb', 'Sftpc', 'Sftpd', 'Lamp3', 'Krt8', 'Cdkn1a', 'Hopx', 'Rtkn2', 'Ager', 'Sox8' , 'Calca', 'Mki67')
DotPlot(fig1, features = epi_genes_figure, cluster.idents= TRUE, dot.scale = 6) +
  RotatedAxis() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

DotPlot(fig1, features = c('Hif1a', 'Epas1'), cluster.idents= T, dot.scale = 6) +
  RotatedAxis() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


library(alluvial)


alluvial_ts(fig1_alluvial,  wave = 15, ygap = 5,
            plotdir = 'centred', alpha=.9,
            grid = TRUE, grid.lwd = 0.2, xmargin = 0.5, lab.cex = 1, xlab = '',
            ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
            leg.col='white', 
            title = "Time Course Cell Populations ")

#MSigDB Hypoxia gene list imported
fig1 <- AddModuleScore(
  object = fig1,
  features = list(mus_hypoxia$gene),
  ctrl = 100,
  name = 'HypoxiaaScore',
  seed = 1234
)



Scores <- list('HIFscore' = fig1@meta.data$HypoxiaaScore1, 'CellCluster' = fig1@meta.data$celltype2, 'Group' = fig1@meta.data$orig.ident )
write.csv(Scores, file = 'Mouse_HIF_Module_scoringXcluster.csv')



pairwise.wilcox.test(x = Scores$HIFscore, g =  Scores$CellCluster, data = Scores, conf.int =T)
pairwise.wilcox.test(x = Scores$HIFscore, g =  Scores$Group, data = Scores, conf.int =T)
kruskal.test(Scores$HIFscore ~ Scores$Group, data = Scores)
wilcox.test(HIFscore ~ Group, data = Scores, conf.int =T, estimate = T, subset = Group %in% c("Unchallenged", "Rep Bleo"))
wilcox.test(HIFscore ~ Group, data = Scores, conf.int =T, estimate = T, subset = Group %in% c("Unchallenged", "SD Bleo"))
pairwise.wilcox.test(HIFscore ~ CellCluster, data = Scores, conf.int =T, estimate = T, subset = Group %in% c("Rep Bleo"))

#ChIP-seq targets from both HIF forms converted then imported. Derived from same set as Alveolar organoids (2011 Bono paper)
adata <- AddModuleScore(
  object = adata,
  features = list(Epas1_gene_targets$symbol),
  ctrl = 100,
  name = 'Epas1Score',
  seed = 1234
)
adata <- AddModuleScore(
  object = adata,
  features = list(Hif1a_gene_targets$symbol),
  ctrl = 100,
  name = 'Hif1aScore',
  seed = 1234
)

Scores_total <- list('GenericHypoxiaScore' = adata@meta.data$HypoxiaaScore1, 'HIF1ascore' = adata@meta.data$Hif1aScore1,'Epas11ascore' = adata@meta.data$Epas1Score1 ,'CellCluster' = adata@meta.data$celltype2, 'Group' = adata@meta.data$orig.ident )
write.csv(Scores_total, file = 'july_2024_Mouse_completeHIFmodule_scoringXcluster.csv')



###Random plot attempts for supplemental figures focusing on hypoxia scoring over time
Idents(fig1) <- 'orig.ident'
fig1_sub <- subset(fig1, idents = c('Unchallenged','Rep Bleo'))
Idents(fig1_sub) <- 'orig.ident'
#fixing the ordering of the plot
Idents(fig1_sub, cells = WhichCells(fig1_sub, idents = c("Rep Bleo"))) <- "Rep Bleo"
Idents(fig1_sub, cells = WhichCells(fig1_sub, idents = c("Unchallenged"))) <- "Unchallenged"
fig1_sub$reorder <- Idents(fig1_sub)
Idents(fig1_sub) <- 'reorder'

fig1_sub <- RunHarmony(fig1_sub, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
fig1_sub <- FindNeighbors(fig1_sub, reduction = 'harmony', dims = 1:45)
fig1_sub <- RunUMAP(fig1_sub, dims = 1:45, reduction = 'harmony')


Idents(fig1_sub) <- 'celltype2'
VlnPlot(fig1_sub, features = c("HypoxiaaScore1"), 
        idents = c('AT1', 'AT2','Alveolar intermediate','BASC-like', 'MCC' ,'Proliferating'), 
        pt.size = 0.5, 
        split.by = 'reorder',
        split.plot = T,
        )

VlnPlot(fig1_sub, features = c("HypoxiaaScore1"), 
        pt.size = 1, 
        split.by = 'reorder',
        split.plot = T,
)

###SD and rep bleo call out plots for response to reviewer
Idents(adata) <- 'orig.ident'
fig1_letter <- subset(adata, idents = c('Unchallenged', 'SD Bleo' ,'Rep Bleo'))
Idents(fig1_letter) <- 'orig.ident'
#fixing the ordering of the plot
Idents(fig1_letter, cells = WhichCells(fig1_letter, idents = c("Rep Bleo"))) <- "Rep Bleo"
Idents(fig1_letter, cells = WhichCells(fig1_letter, idents = c("SD Bleo"))) <- "SD Bleo"
Idents(fig1_letter, cells = WhichCells(fig1_letter, idents = c("Unchallenged"))) <- "Unchallenged"
fig1_letter$reorder <- Idents(fig1_letter)
Idents(fig1_letter) <- 'reorder'

Idents(fig1_letter) <- 'celltype2'
VlnPlot(fig1_letter, features = c("HypoxiaaScore1"), 
        pt.size = 0.1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF")),
        add.noise = F
        )

VlnPlot(fig1_letter, features = c("Hif1aScore1"), 
        pt.size = 1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF"))
)

VlnPlot(fig1_letter, features = c("Epas1Score1"), 
        pt.size = 1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF"))
)

VlnPlot(fig1_letter, features = c("HypoxiaaScore1"), 
        idents = c('AT1', 'AT2','Alveolar intermediate','BASC-like', 'MCC' ,'Proliferating'),
        pt.size = 0.1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF"))
)

VlnPlot(fig1_letter, features = c("Hif1aScore1"), 
        idents = c('AT1', 'AT2','Alveolar intermediate','BASC-like', 'MCC' ,'Proliferating'),
        pt.size = 1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF"))
)

VlnPlot(fig1_letter, features = c("Epas1Score1"), 
        idents = c('AT1', 'AT2','Alveolar intermediate','BASC-like', 'MCC' ,'Proliferating'),
        pt.size = 1, 
        split.by = 'reorder',
        cols = c(c("#757576", "#C49A6C", "#FFFFFF"))
)

##Other Subplots 

VlnPlot(fig1_sub, features = c("Hif1aScore1", "Epas1Score1"), 
        idents = c('AT1', 'AT2','Alveolar intermediate','BASC-like', 'MCC' ,'Proliferating'), 
        pt.size = 0.5, 
        split.by = 'reorder',
        split.plot = T,
)

table(Idents(fig1_sub),fig1_sub$celltype2)

DimPlot_scCustom(fig1_sub,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6.0, colors_use = c("#757576", "#010101"))
Idents(fig1_sub) <- 'celltype2'
DimPlot_scCustom(fig1_sub,  label = F, figure_plot = TRUE , raster = TRUE, pt.size = 6.0, colors_use = col_var)
Plot_Density_Custom(seurat_object = fig1_sub, features = c("Hif1aScore1"))
Plot_Density_Custom(seurat_object = fig1_sub, features = c("Epas1Score1"))

Idents(fig1_sub) <- 'celltype2'
VlnPlot_scCustom(fig1_sub, features = "HypoxiaaScore1", colors_use = col_var) & NoLegend()


pairwise.wilcox.test(x = Scores$HIFscore, g =  Scores$CellCluster, data = Scores, conf.int =T)
pairwise.wilcox.test(x = Scores$HIFscore, g =  Scores$Group, data = Scores, conf.int =T)
kruskal.test(Scores$HIFscore ~ Scores$Group, data = Scores)
wilcox.test(HIFscore ~ Group, data = Scores, conf.int =T, estimate = T, subset = Group %in% c("Unchallenged", "Rep Bleo"))
wilcox.test(HIFscore ~ CellCluster, data = Scores, conf.int =T, estimate = T, subset = CellCluster %in% c("AT2", "AT1"))

#HIF KO workup for Figures____________________________________________________________________________________________________________________________________________________________________
table(Idents(fig2),fig2$orig.ident)

Idents(fig2) <- 'orig.ident'
my_levels <- c('Unchallenged', 'Rep Bleo', 'Hif1/2 Unchallenged', 'Hif1/2 Rep Bleo')
Idents(fig2) <- factor(Idents(fig2), levels = my_levels)

Idents(fig2) <- 'orig.ident'
DimPlot_scCustom(fig2, figure_plot = T, raster = TRUE, pt.size = 6, label.size = 7, label = F, colors_use = c("#757575","#FB0106","#000000","#0000FF"))

Idents(fig2) <- 'celltype2'
DimPlot_scCustom(fig2,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))


#loading prop table formatted 
#fig2_WT_allucial <- read_excel("revamp_prop_table_20230718_for_alluvial.xlsx", 
#+     sheet = "Fig2 WT")
#+     
# fig2_KO_allivial <- read_excel("revamp_prop_table_20230718_for_alluvial.xlsx", 
#+     sheet = "Fig2-HIF KO")

WTrep<- alluvial_ts(fig2_WT_allucial,  wave = 15, ygap = 5,
            plotdir = 'centred', alpha=.9,
            grid = TRUE, grid.lwd = 0.2, xmargin = 0.5, lab.cex = 1, xlab = '',
            ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
            leg.col='white', 
            title = "WT Rep Bleo")

KOrep <- alluvial_ts(fig2_KO_allivial,  wave = 15, ygap = 5,
                     plotdir = 'centred', alpha=.9,
                     grid = TRUE, grid.lwd = 0.2, xmargin = 0.5, lab.cex = 1, xlab = '',
                     ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
                     leg.col='white', 
                     title = "KO Rep Bleo")

#Cemitool analysis
library(CEMiTool)

##https://www.bioconductor.org/packages/devel/bioc/vignettes/CEMiTool/inst/doc/CEMiTool.html#gene-filtering

#HIF KO workup for Figures- Reassigning genotypes for sets
Idents(fig2) <- 'orig.ident'
Idents(fig2, cells = WhichCells(fig2, idents = c('Unchallenged', 'Rep Bleo'))) <- "WT"
Idents(fig2, cells = WhichCells(fig2, idents = c('Hif1/2 Unchallenged', 'Hif1/2 Rep Bleo'))) <- "HIF KO"
fig2$genotype <- Idents(fig2)
Idents(fig2) <- 'genotype'

table(Idents(fig2),fig2$orig.ident)



#subset to only secretory lineages
Idents(fig2) <- 'celltype2'
adata3_sec <- subset(fig2, idents = c('Basal-like', 'Secretory', 'Secretory-Muc5b+', 'BASC-like', 'Airway intermediate','M-like cells'))
Idents(adata3_sec) <- 'genotype'
table(Idents(adata3_sec),adata3_sec$celltype2)

counts.ko.sec <- as.data.frame(adata3_sec@assays$SCT@counts)
counts.ko.sec.RNA <- as.data.frame(adata3_sec@assays$RNA@data)
annot.ko.sec <- as.data.frame(adata3_sec@active.ident)
Idents(adata3_sec) <- 'celltype2'
annot.ko.sec_celltype <- as.data.frame(adata3_sec@active.ident)

#writeeCSV so that colulnm 1 can be named "SampleName" and column with genotype can be named "Class", then reimport
write.csv(annot.ko.sec, file='sec_annot_for_KO.csv')
write.csv(annot.ko.sec_celltype, file='sec_annot_for_KO_celltype2.csv')

adata3_sec <- RunPCA(adata3_sec, verbose = F)
adata3_sec <- RunHarmony(adata3_sec, group.by.vars = 'orig.ident', max.iter.harmony = 50, theta = 1, nclust = 25)
adata3_sec <- FindNeighbors(adata3_sec, reduction = 'harmony', dims = 1:45)
adata3_sec <- RunUMAP(adata3_sec, dims = 1:45, reduction = 'harmony')
Idents(adata3_sec) <- 'celltype2'
sec_colors <- c("#80FF00","#5CCC5C", "#454599", "#A64DFF", "#CC1FCC", "#99004D" )
DimPlot_scCustom(adata3_sec, figure_plot = T, split_seurat = T, pt.size = 2, split.by = 'genotype', label = F,colors_use = sec_colors)


#reimported after as sec_annot_for_KO
                     
reactome_gmt <- read_gmt("C:/Users/mccalas1/OneDrive - VUMC/Desktop/HIF scRNA/CemiTool HIF work/Mouse GMT files/reactome.v5.2.symbols_mouse.gmt")
                     
#import PPI sheet after conversion to mouse gene names and FORMAT AS A DATA FRAME
                     
interactions <- as.data.frame(interactions_filtered)

#apply_vst is needed for scRNA data
cem14<- cemitool(counts.ko.sec, sec_annot_for_KO, gmt = reactome_gmt, interactions = interactions , set_beta = 14, apply_vst = T, filter_pval = 0.15, cor_method = 'spearman', cor_function = 'bicor', gsea_max_size = 1500,  plot = T, plot_diagnostics = T, verbose = T)
show_plot(cem14, "gsea")


cem14_celltype<- cemitool(counts.ko.sec, sec_annot_for_KO_celltype2, gmt = reactome_gmt, interactions = interactions , set_beta = 14,apply_vst = T, filter_pval = 0.15, cor_method = 'spearman', cor_function = 'bicor', gsea_max_size = 1500,  plot = T, plot_diagnostics = T, verbose = T)
show_plot(cem14_celltype, "gsea")

cem5<- cemitool(counts.ko.sec, sec_annot_for_KO, gmt = reactome_gmt, interactions = interactions , set_beta = 5, apply_vst = T, filter_pval = 0.15, cor_method = 'spearman', cor_function = 'bicor', gsea_max_size = 1500,  plot = T, plot_diagnostics = T, verbose = T)
show_plot(cem5, "gsea")

cem5_celltype<- cemitool(counts.ko.sec, sec_annot_for_KO_celltype2, gmt = reactome_gmt, interactions = interactions , set_beta = 5,apply_vst = T, filter_pval = 0.15, cor_method = 'spearman', cor_function = 'bicor', gsea_max_size = 1500,  plot = T, plot_diagnostics = T, verbose = T)
show_plot(cem5_celltype, "gsea")

cem5 <- mod_gsea(cem = cem5, gsea_max_size = 10000, gsea_min_size=2, verbose = T)

library(ggplot2)
options(ggrepel.max.overlaps = Inf)
cem5 <- plot_interactions(cem5, n=15) # generate plot
plots <- show_plot(cem5, "interaction") #
                     plots[1]
                     plots[2]
                     plots[3]
                     plots[4]
                     plots[5]
                     plots[6]
                     
hubs <- get_hubs(cem5,50)
save_plots(cem5, value="interaction", directory="./Plots_Cemitool_July_v4_beta5_15PPI_plot")
                     
# write analysis results into files
write_files(cem14, directory="./Tables_CemiTool_July_v4")
                     
# save all plots
save_plots(cem14, "all", directory="./Plots_CemiTool_July_v4")

# create report as html document
generate_report(cem14, directory="./Report_CemiTool_July_v4")

# create report as html document
generate_report(cem14_celltype, directory="./Report_CemiTool_July_v4_celltype")

# write analysis results into files
write_files(cem14_celltype, directory="./Tables_CemiTool_July_v4_celltype")

# save all plots
save_plots(cem14_celltype, "all", directory="./Plots_CemiTool_July_v4_celltype")
                     
                     #writing top 30 genes for hubs
                     write.csv(hubs$M1, file = 'Hubs1_for_final_modules.csv')
                     write.csv(hubs$M2, file = 'Hubs2_for_final_modules.csv')
                     write.csv(hubs$M3, file = 'Hubs3_for_final_modules.csv')
                     write.csv(hubs$M4, file = 'Hubs4_for_final_modules.csv')
                     write.csv(hubs$M5, file = 'Hubs5_for_final_modules.csv')
                     write.csv(hubs$M6, file = 'Hubs6_for_final_modules.csv')
                     write.csv(hubs$M7, file = 'Hubs7_for_final_modules.csv')
                     write.csv(hubs$M8, file = 'Hubs8_for_final_modules.csv')
                     write.csv(hubs$Not.Correlated, file = 'Hubs_nonCorr_for_final_modules.csv')
                     
                     summary <- mod_summary(cem14, method = c('mean','median', 'eigengene'), verbose = TRUE)
                     
                     
                     # write analysis results into files
                     write_files(cem5, directory="./Tables_CemiTool_July_v4_beta5")
                     
                     # save all plots
                     save_plots(cem5, "all", directory="./Plots_CemiTool_July_v4_beta5")
                     
                     # create report as html document
                     generate_report(cem5, directory="./Report_CemiTool_July_v4_beta5")
                     
                     # create report as html document
                     generate_report(cem5_celltype, directory="./Report_CemiTool_July_v4_celltype_beta5")
                     
                     # write analysis results into files
                     write_files(cem5_celltype, directory="./Tables_CemiTool_July_v4_celltype_beta5")
                     
                     # save all plots
                     save_plots(cem5_celltype, "all", directory="./Plots_CemiTool_July_v4_celltype_beta5")
                     
                     #writing top 30 genes for hubs
                     write.csv(hubs$M1, file = 'Hubs1_for_final_modules.csv')
                     write.csv(hubs$M2, file = 'Hubs2_for_final_modules.csv')
                     write.csv(hubs$M3, file = 'Hubs3_for_final_modules.csv')
                     write.csv(hubs$M4, file = 'Hubs4_for_final_modules.csv')
                     write.csv(hubs$M5, file = 'Hubs5_for_final_modules.csv')
                     write.csv(hubs$M6, file = 'Hubs6_for_final_modules.csv')
                     write.csv(hubs$M7, file = 'Hubs7_for_final_modules.csv')
                     write.csv(hubs$M8, file = 'Hubs8_for_final_modules.csv')
                     write.csv(hubs$Not.Correlated, file = 'Hubs_nonCorr_for_final_modules.csv')
                     
                     summary <- mod_summary(cem14, method = c('mean','median', 'eigengene'), verbose = TRUE)

saveRDS(cem5, file = "Cemitool_Beta5_secretoryFinal.rds")
saveRDS(cem5_celltype, file = "Cemitool_Beta5_secretoryFinal_celltype.rds")
saveRDS(cem14, file = "Cemitool_Beta14_secretoryFinal.rds")
saveRDS(cem14_celltype, file = "Cemitool_Beta14_secretoryFinal_celltype.rds")


                     

#enrichR plots for Cemitool moduls
                     
modules_gene <- read_excel("Tables_CemiTool_July_v4_beta5/modules_genes_pathway_formatted.xlsx", sheet = "Sheet1")
dbs <- listEnrichrDbs()
qdbs <- c("MSigDB_Hallmark_2020", "ARCHS4_TFs_Coexp", "PPI_Hub_Proteins", "KEA_2015")

enriched2 <- enrichr(c(modules_gene$M2), qdbs)
enriched3 <- enrichr(c(modules_gene$M3), qdbs)
enriched4 <- enrichr(c(modules_gene$M4), qdbs)
enriched5 <- enrichr(c(modules_gene$M5), qdbs)
enriched6 <- enrichr(c(modules_gene$M6), qdbs)
enriched8 <- enrichr(c(modules_gene$M8), qdbs)

plotEnrich(enriched2[[1]],showTerms = 10, y = "Ratio")
plotEnrich(enriched3[[1]],showTerms = 10, y = "Ratio")
plotEnrich(enriched4[[1]],showTerms = 10, y = "Ratio")
plotEnrich(enriched5[[1]],showTerms = 10, y = "Ratio")
plotEnrich(enriched6[[1]],showTerms = 10, y = "Ratio")
plotEnrich(enriched8[[1]],showTerms = 10, y = "Ratio")

plotEnrich(enriched2[[1]],showTerms = 5, y = "Ratio")
plotEnrich(enriched4[[1]],showTerms = 5, y = "Ratio")
plotEnrich(enriched5[[1]],showTerms = 5, y = "Ratio")
plotEnrich(enriched6[[1]],showTerms = 5, y = "Ratio")
plotEnrich(enriched8[[1]],showTerms = 5, y = "Ratio")

plotEnrich(enriched2[[1]],showTerms = 10, y = "p-value")
plotEnrich(enriched4[[1]],showTerms = 10, y = "p-value")
plotEnrich(enriched5[[1]],showTerms = 10, y = "p-value")
plotEnrich(enriched6[[1]],showTerms = 10, y = "p-value")
plotEnrich(enriched8[[1]],showTerms = 10, y = "p-value")

plotEnrich(enriched2[[1]],showTerms = 5, y = "p-value")
plotEnrich(enriched4[[1]],showTerms = 5, y = "p-value")
plotEnrich(enriched5[[1]],showTerms = 5, y = "p-value")
plotEnrich(enriched6[[1]],showTerms = 5, y = "p-value")
plotEnrich(enriched8[[1]],showTerms = 5, y = "p-value")

plotEnrich(enriched2[[3]],showTerms = 5, y = "Ratio")
plotEnrich(enriched4[[3]],showTerms = 5, y = "Ratio")
plotEnrich(enriched5[[3]],showTerms = 5, y = "Ratio")
plotEnrich(enriched6[[3]],showTerms = 5, y = "Ratio")
plotEnrich(enriched8[[3]],showTerms = 5, y = "Ratio")

plotEnrich(enriched2[[3]],showTerms = 10, y = "p-value")
plotEnrich(enriched4[[3]],showTerms = 10, y = "p-value")
plotEnrich(enriched5[[3]],showTerms = 10, y = "p-value")
plotEnrich(enriched6[[3]],showTerms = 10, y = "p-value")
plotEnrich(enriched8[[3]],showTerms = 10, y = "p-value")

plotEnrich(enriched2[[4]],showTerms = 10, y = "p-value")
plotEnrich(enriched4[[4]],showTerms = 10, y = "p-value")
plotEnrich(enriched5[[4]],showTerms = 10, y = "p-value")
plotEnrich(enriched6[[4]],showTerms = 10, y = "p-value")
plotEnrich(enriched8[[4]],showTerms = 10, y = "p-value")

#_figure 3 LINEAGE TRACE FIGURE___________________________________________________________________________________________________________________________________________________________________________________________________________
table(Idents(fig3),fig3$orig.ident)

Idents(fig3) <- 'orig.ident'
my_levels_lin <- c('Scgb1a1-trace', 'Rep Bleo', 'Sftpc-trace')
Idents(fig3) <- factor(Idents(fig3), levels = my_levels_lin)

Idents(fig3) <- 'orig.ident'
DimPlot_scCustom(fig3, figure_plot = T, raster = TRUE, pt.size = 6, label.size = 7, label = F, colors_use = c("#92298D", "#000000", "#369846"))

Idents(fig3) <- 'celltype2'
DimPlot_scCustom(fig3,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))
DimPlot_scCustom(fig2,  label = F, figure_plot = TRUE, raster = TRUE, label.size = 5,  pt.size = 1.0, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))
DimPlot_scCustom(fig2,  label = T, figure_plot = TRUE, raster = TRUE, label.size = 5,  pt.size = 1.0, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow")),"#000000","#0000FF"))

Idents(fig2) <- 'celltype2'
DimPlot_scCustom(fig2,  label = F, figure_plot = TRUE, raster = TRUE, pt.size = 6, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))
DimPlot_scCustom(fig2,  label = F, figure_plot = TRUE, raster = TRUE, label.size = 5,  pt.size = 1.0, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))
DimPlot_scCustom(fig2,  label = T, figure_plot = TRUE, raster = TRUE, label.size = 5,  pt.size = 1.0, colors_use = DiscretePalette_scCustomize(num_colors = 12, palette = "varibow"))

LINEAGE_REP <- alluvial_ts(revamp_prop_table_20230718,  wave = 15, ygap = 5,
                     plotdir = 'centred', alpha=.9,
                     grid = TRUE, grid.lwd = 0.2, xmargin = 0.5, lab.cex = 1, xlab = '',
                     ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
                     leg.col='white', 
                     title = "Scgb1a1 vs Reo vs Sftpc")







                 

                 
                 
                 