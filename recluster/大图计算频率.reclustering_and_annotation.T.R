library(Seurat, lib.loc = "/home/data/refdir/Rlib")
seu=readRDS("~/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/RPCA.CUIMC/RPCA.CUIMC_int_obj.ScaleData.RunPCA.FindClusters.rds")
Idents(seu) <- seu$integrated_snn_res.0.1
data.seurat.filtered <- RenameIdents(seu,
                                     `0` = "T/NK",
                                     `1` = "Epithelia",
                                     `2` = "Fibroblasts",
                                     `3` = "Myeloid cells",
                                     `4` = "B/plasma",
                                     `5` = "Endothelium",
                                     `6` = "B/plasma",
                                     `7` = "Epithelia",
                                     `8` = "Pericytes",
                                     `9` = "pDCs")
consistentcolors=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#34bfd6","#D65C06","#CC78AA","#56B4EA","#E59C04","#0868B6","#94A1A1","#02996E","#EEE53F")
consistentcolors2 = colors <- c("#006E82", "#AA0A3C")
DimPlot(data.seurat.filtered,reduction = "umap",cols = consistentcolors,raster=FALSE)
rm(seu)
DimPlot(scRNA_harmony,reduction = "umap",cols = consistentcolors,raster=FALSE)
# c("#F0C6BE","#DFD57C","#DEC9BC","#A6C694","#5DB7D2","#42A6E5","#CDE0E8","#506AD1","#D1AAEF","#2F8EA1",))
# DimPlot(scRNA_harmony,reduction = "tsne",split.by = "timepoint",cols = c("#F0C6BE","#DFD57C","#DEC9BC","#A6C694","#5DB7D2","#42A6E5","#CDE0E8","#506AD1","#D1AAEF","#2F8EA1"))
pdf(file="/home/data/t180330/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/Major clusters and respective cell-type assignmentsplit.bygroup.umap.pdf", width = 12, height = 5)
DimPlot(scRNA_harmony,split.by = "group",reduction = "umap",cols = consistentcolors,raster=FALSE,label=F)
dev.off()
pdf(file="/home/data/t180330/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/Major clusters and respective cell-type assignment.umap.pdf", width = 8, height = 5)
DimPlot(scRNA_harmony,reduction = "umap",cols = consistentcolors,raster=FALSE,label=F)
dev.off()
pdf(file="/home/data/t180330/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/Major clusters and respective cell-type assignmentsplit.bygroupwith label.umap.pdf", width = 12, height = 5)
DimPlot(scRNA_harmony,split.by = "group",reduction = "umap",cols = consistentcolors,raster=FALSE,label=T)
dev.off()
pdf(file="/home/data/t180330/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/Major clusters and respective cell-type assignment.umapwith label.pdf", width = 8, height = 5)
DimPlot(scRNA_harmony,reduction = "umap",cols = consistentcolors,raster=FALSE,label=T)
dev.off()
pdf(file="/home/data/t180330/202311_tnbc_caf_05_20240521/DATA/MERGE.ALL/Major clusters and respective cell-type assignment.umapsplit.group.pdf", width = 6.5, height = 5)
DimPlot(scRNA_harmony,reduction = "umap",group.by = "group",cols = consistentcolors2,raster=FALSE,label=F)
dev.off()
table(scRNA_harmony@meta.data$group)
write.csv(table(scRNA_harmony@active.ident), '~/202311_tnbc_caf_05_20240521/DATA//table_freq_cell_type_main.csv', row.names = F)
write.csv(prop.table(table(subset(scRNA_harmony, group == 'Tumor')@active.ident)), '~/202311_tnbc_caf_05_20240521/DATA//table_freq_cell_type_main_tummor.csv', row.names = F)
write.csv(prop.table(table(subset(scRNA_harmony, group == 'Normal')@active.ident)), '~/202311_tnbc_caf_05_20240521/DATA/table_freq_cell_type_main_ctr.csv', row.names = F)

#
sel.genes = c("KRT19","KRT17","KRT18","COL1A1","LUM","DCN",'CD3D',"CD3E","LYZ","CD68","VWF","PECAM1","CD79A","CD38","IGKC","ACTA2","MCAM","IRF7","GZMB")

p <- DotPlot(object = scRNA_harmony, features = sel.genes, assay = "RNA", scale = T) + 
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  scale_x_discrete(breaks = sel.genes, labels = sel.genes)
p


data_lungs_all=scRNA_harmony
scRNA_harmony@active.ident
table(data_lungs_all$TIG)
df_tobesummed = data.frame(orig.ident = data_lungs_all$TIG, group = data_lungs_all$group, 
                           cell_type_fine = data_lungs_all@active.ident, immune_status = data_lungs_all$immune_status)
df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_fine, 
                                       immune_status, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
setwd("~/202311_tnbc_caf_05_20240521/DATA")
write.table(df_summed, "boxplot_proportions_fine.csv", sep = ",", row.names = F, 
            col.names = T, quote = F)
# calculate frequencies of cell_type_main classes in each sample
data_lungs_all <- seu_filtered2

df_tobesummed_main = data.frame(orig.ident = data_lungs_all$TIG, group = data_lungs_all$group, 
                                cell_type_fine = data_lungs_all@active.ident, immune_status = data_lungs_all$immune_status)
df_summed_main = df_tobesummed_main %>% group_by(orig.ident, cell_type_fine, immune_status, 
                                                 group) %>% tally()
df_summed_main = df_summed_main %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
ggboxplot(df_summed_main, x = "cell_type_fine", y = "freq", color = "group", add = "jitter") + 
  ylim(0, 0.8) + stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors2)
ggsave(paste0("Figure_1E.pdf"), width = 11, height = 7)

saveRDS(data_lungs_all,"~/202311_tnbc_caf_05_20240521/DATA/rpca.all.main.ident.rds")


###!!!~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/reclustering_and_annotation.T.R
#### Fine cell type annotation of B cell subset for MBM_sc
#### Author: Jana Biermann, PhD
# load("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/20230823.rds.RData")
lib_loc <-  "/home/data/refdir/Rlib"
library(Seurat,lib.loc = lib_loc)
library(devtools,lib.loc = lib_loc)
# install_github("YosefLab/VISION")
library(VISION)
# devtools::install_github("caleblareau/BuenColors",force = TRUE)
library(BuenColors)
set.seed(1234)

signatures <- "/home/data/t180330/h.all.v6.2.symbols.gmt"
vision.obj <- Vision(data.seurat.filtered, signatures = signatures, pool = F)
# options(mc.cores = 10)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)

data.seurat.filtered[['GSEA']] <- CreateAssayObject(counts = t(sigScores))
DefaultAssay(data.seurat.filtered) <- "GSEA"



saveRDS(data.seurat.filtered,"data.seurat.filtered.GSEA.RDS")



library(dplyr)
lib_loc="/home/data/refdir/Rlib"
library(Seurat,lib.loc=lib_loc)
library(ggplot2)
library(gplots)
# library(viridis,lib.loc=lib_loc)
'%notin%' <- Negate('%in%')
path.ct <- '/home/data/t180330/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts'
filename <- 'Fibroblasts'
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
library(future)
plan("multicore", workers = 6)
options(future.globals.maxSize = 4000 * 1024^6)
seu <- subset(data.seurat.filtered, seurat_clusters == 'Fibroblasts')
seu <- ScaleData(object = seu)
seu <- FindVariableFeatures(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(object = seu, dims = 1:30)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'ident', shuffle = T, raster = T)+DimPlot(seu, reduction = 'umap', label = T, group.by = 'Tissue', shuffle = T, raster = T)
#####重新整合########

###重新整合CCAF###
seu <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/seu.withcluster-rename.rds")
seu$seurat_clusters <- Idents(seu)
DefaultAssay(seu) <- 'RNA'
CAF=subset(seu,seurat_clusters=="T/NK")
sf.list <- SplitObject(seu,split.by="TIG")
sf.list <- lapply(X = sf.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sf.list,nfeatures = 1000) #Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance
  s.genes=Seurat::cc.genes.updated.2019$s.genes
  g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
features <- setdiff(features,c(s.genes,g2m.genes))  ##  Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance
# install.packages("future_eapply")  future_lapply
sf.list <- lapply(X = sf.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = sf.list, anchor.features = features, reduction = "rpca")
nCoV.integrated <- IntegrateData(anchorset = immune.anchors)

saveRDS(nCoV.integrated,"fib.integrated.rds")
 saveRDS(seu,"fib.RDS")
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, "markers.csv", row.names = F)
DefaultAssay(seu) <- "RNA"
FeaturePlot(seu,features = c("CD4", "CD8A", "FOXP3", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "GZMB", "PRF1", "SELL", "KLRG1", "KLF2"))
DimPlot(object = seu, reduction = 'umap',label = TRUE,label.size = 8,repel = T,cols=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#F69B9A",
                                                                                                       "#D65C06","#CC78AA"))
DimPlot(object = seu, reduction = 'umap',label = TRUE,label.size = 8,repel = T,split.by = "tissue",cols=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#F69B9A",
                                                                                                                          "#D65C06","#CC78AA"))

saveRDS(seu,"seu.wtih.cluster.rds") 
# data.seurat.filtered <- RenameIdents(seu,
#                                      `0` = "CD4-CCR7",
#                                      `1` = "CD8-CXCL13-GZMK",
#                                      `2` = "CD4-CXCL13",
#                                      `3` = "Treg-FOXP3",
#                                      `4` = "CD8-IL7R",
#                                      `5` = "NK",
#                                      `6` = "Proliferating T",
#                                      `7` = "CD8-γδT-TRDC",
#                                      `8` = "Epithelia"
data.seurat.filtered <- RenameIdents(seu,
                                     `0` = "CD4",
                                     `1` = "CD8-CXCL13",
                                     `2` = "CD4",
                                     `3` = "Treg-FOXP3",
                                     `4` = "CD8-KLRC1",
                                     `5` = "NK",
                                     `6` = "Proliferating T",
                                     `7` = "CD8-γδT-TRDC",
                                     `8` = "Epithelia")
data.seurat.filtered <- subset(data.seurat.filtered,cell_type_fine!="Epithelia")
CD8 <- subset(data.seurat.filtered,cell_type_fine %in% c("CD8-CXCL13","CD8-KLRC1","CD8-γδT-TRDC"))
# A
DimPlot(CD8, label = T)

# B
gs <- c("GPM6B", "IGFBP7", "DAPL1", "BDH2", "CXCR3", "ITGB7", "CCR2", "GZMK", "GZMF", "GZMC",
        "FCER1G", "GZMB", "LAG3", "HAVCR2", "BHLHE40", "PRF1", "MCM3", "MCM5", "MCM6", "LIG1", "STMN1",
        "BIRC5", "PCLAF", "MKI67", "CCL3", "CCL4", "IFNG", "CCL1", "IFIT3", "IFIT1", "CXCL10", "ISG15")
avg <- AverageExpression(CD8, "RNA", gs)
ep <- t(scale(t(avg$RNA)))
pheatmap::pheatmap(ep, cluster_rows = F, cluster_cols = F)

markers <- FindAllMarkers(data.seurat.filtered, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers,'markers_fine_celltype.csv', row.names = F)
data.seurat.filtered$cell_type_fine <- Idents(data.seurat.filtered)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops
DefaultAssay(data.seurat.filtered) <- 'RNA'
data.seurat.filtered <- NormalizeData(data.seurat.filtered) %>%
  ScaleData()
DefaultAssay(data.seurat.filtered) <- 'integrated'

pdf(file = 'markers_heatmap_celltype.pdf', width = 10,
    height = 18)
DoHeatmap(data.seurat.filtered, features = tops$gene, group.by = 'cell_type_fine', raster = T,
          assay = 'RNA')
dev.off()



data_lungs_all=seu
table(data_lungs_all$TIG)
df_tobesummed = data.frame(orig.ident = data_lungs_all$TIG, group = data_lungs_all$group, 
                           cell_type_fine = data_lungs_all@active.ident, immune_status = data_lungs_all$immune_status)
df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_fine, 
                                       immune_status, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
setwd("~/202311_tnbc_caf_05_20240521/DATA")
write.table(df_summed, "boxplot_proportions_fine.csv", sep = ",", row.names = F, 
            col.names = T, quote = F)
# calculate frequencies of cell_type_main classes in each sample
df_tobesummed_main = data.frame(orig.ident = data_lungs_all$TIG, group = data_lungs_all$group, 
                                cell_type_fine = data_lungs_all@active.ident, immune_status = data_lungs_all$immune_status)
df_summed_main = df_tobesummed_main %>% group_by(orig.ident, cell_type_fine, immune_status, 
                                                 group) %>% tally()
df_summed_main = df_summed_main %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
ggboxplot(df_summed_main, x = "cell_type_fine", y = "freq", color = "group", add = "jitter") + 
  ylim(0, 0.8) + stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors2)
ggsave(paste0("Figure_1E.pdf"), width = 11, height = 7)

saveRDS(data_lungs_all,"~/202311_tnbc_caf_05_20240521/DATA/rpca.all.main.ident.rds")



# for (i in 1:length(sf.list)) {
#   sf.list[[i]] <- NormalizeData(sf.list[[i]], verbose = FALSE)
#   sf.list[[i]] <- FindVariableFeatures(sf.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
# 
# 
# seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes,
#                         g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
# DefaultAssay(seu) <- 'integrated'
# seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling',
#                          'non-cycling')
# # Differential gene expression (DGE) based on clusters
# markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
#                           logfc.threshold = 0.25, test.use = 'MAST')
# write.csv(markers, paste0(path.ct, 'markers_', filename, '.csv'), row.names = F)
