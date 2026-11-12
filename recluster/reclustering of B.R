#### Fine cell type annotation of B cell subset for MBM_sc
#### Author: Jana Biermann, PhD
data.seurat.filtered <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/seu.withcluster-rename.rds")
data.seurat.filtered$cell_type_fine <- Idents(data.seurat.filtered)
# markers <- FindAllMarkers(data.seurat.filtered, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
#                           logfc.threshold = 0.25, test.use = 'MAST')
# write.table(markers,"all.eightcluster.marker.tsv",sep = '\t')
setwd("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B")
library(dplyr)
# lib_loc="/home/data/refdir/Rlib"
library(Seurat,lib.loc="~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/spatial/libPaths")
library(SeuratObject,lib.loc="~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/spatial/libPaths")
library(ggplot2)
library(gplots)
library(viridis)
devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)
install.packages('Rcpp')
detach("package:Rcpp", unload=TRUE)

library(Rcpp)
install.packages('harmony')
install.packages("devtools")
library(devtools)
update_packages()

#使用了Seurat V4
seu <- subset(data.seurat.filtered, cell_type_fine == 'B/plasma')
DefaultAssay(seu) <- 'RNA'
#### Remove IG genes ####
seu <- seu[grep("^IGJ|^IGH|^IGK|^IGL", rownames(seu), value = T, invert = T), ]
sce <- seu
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce=CellCycleScoring(object = sce,
                     s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
sce[['percent.mito']] <- PercentageFeatureSet(sce, pattern = "^MT-")
library(future)
plan("multicore", workers = 6)
options(future.globals.maxSize = 4000 * 1024^6) 
sce <- SCTransform(sce, vars.to.regress =c("percent.mito","S.Score","G2M.Score","nCount_RNA"),verbose = TRUE)

# https://github.com/const-ae/glmGamPoi


harmonized_seurat <- RunHarmony(sce, project.dim = FALSE,
                                group.by.vars = c("TIG","tissue"),
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(harmonized_seurat)

x<-c(1:15)
# eigs <- harmonized_seurat[["pca"]]@stdev^2 #PCA的Standard dev的平方#

print(paste0("Variance captured by ",x," PCs: ",sum(eigs[x] / sum(eigs))))
# 30 PCs: 0.864842918535045"
harmonized_seurat <- FindNeighbors(harmonized_seurat, reduction = "harmony", dims = x)
# 5: closing unused connection 5 (<-localhost:11192) 
# Use the integrated_snn graph for clustering
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = seq(0.1, 1.2, 0.1), graph.name = "integrated_snn")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony",assay = "SCT",dims = x)
clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(harmonized_seurat@meta.data),value = TRUE),
                            function(x) length(unique(harmonized_seurat@meta.data[,x])))
clusters_resolution

Idents(harmonized_seurat) <- "integrated_snn_res.0.5"
library(dplyr)
markers <- FindAllMarkers(harmonized_seurat, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
# 假设 markers 是您的数据框
# 使用 dplyr 包按 group 分组并按 avg_log2FC 降序排序
markers$cluster <- as.integer(markers$cluster)
sorted_markers <- markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC))
library(openxlsx)
# 创建一个新的Excel文件
wb <- createWorkbook()

# 为每个cluster创建一个sheet并写入数据
for (cluster in unique(sorted_markers$cluster)) {
  # 获取当前cluster的数据
  cluster_data <- sorted_markers %>%
    filter(cluster == cluster)
  
  # 添加一个新的sheet并写入数据
  openxlsx::addWorksheet(wb, cluster)
  writeData(wb, sheet = cluster, x = cluster_data)
}

# 保存Excel文件
file_path <- "~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B/cluster_markers.xlsx"
saveWorkbook(wb, file_path, overwrite = TRUE)


write.table(sorted_markrs,"all.sixcluster.res.05.B.marker.tsv",sep = '\t')
saveRDS(harmonized_seurat,"harmonized_seurat.RDS")
Idents(harmonized_seurat)
harmonized_seurat <- RenameIdents("0"="Bn-TCL1A",JCHAIN)
# T <- RenameIdents(T,"0"= "CD4-CCR7", "1"= "CD4-FOXP3", "5"= "CD4-CXCL13", "2")


png(file="6-NKT-tsne.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = harmonized_seurat, reduction = 'tsne',label = TRUE)
dev.off()
dpi = 300
png(file="6-NKT-umap.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = harmonized_seurat, reduction = 'umap',label = TRUE)
dev.off()
saveRDS(seu, file = "6-NKT.rds")
png(file="6-NKT-umap-split-sample.png", width = dpi*20, height = dpi*18, units = "px",res = dpi,type='cairo')
DimPlot(object = harmonized_seurat, reduction = 'umap',label = TRUE, split.by = 'TIG', ncol = 5)
dev.off()

png(file="6-NKT-umap-split-group.png", width = dpi*20, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = harmonized_seurat, reduction = 'umap',label = TRUE, split.by = 'tissue',ncol = 3)
dev.off()

png(file="6-NKT-umap-group-sample.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = harmonized_seurat, reduction = 'umap',label = TRUE, group.by = 'TIG')
dev.off()




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
immune.anchors <- FindIntegrationAnchors(object.list = sf.list, anchor.features = features, reduction = "rpca")
nCoV.integrated <- IntegrateData(anchorset = immune.anchors,k.weight=10)
saveRDS(nCoV.integrated,"b.integrated.rds")
seu <- nCoV.integrated
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)


###first generate data and scaledata in RNA assay
DefaultAssay(seu) <- "RNA"
seu[['percent.mito']] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(seu) <- "integrated"
dpi = 300
png(file="6-qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="6-umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
seu <- RunPCA(seu, verbose = FALSE)
#visulaization pca result
seu <- ProjectDim(object = seu)
png(file="6-pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = seu,ndims = 50)
dev.off()

###cluster
seu <- FindNeighbors(object = seu, dims = 1:20)
seu <- FindClusters(object = seu, resolution = 0.4)

# DimPlot(seu,group.by = "integrated_snn_res.0.4",reduction = )



###tsne
seu <- RunTSNE(object = seu, dims = 1:50)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
png(file="6-NKT-tsne.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = seu, reduction = 'tsne',label = TRUE)
dev.off()
png(file="6-NKT-umap.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = seu, reduction = 'umap',label = TRUE)
dev.off()
saveRDS(seu, file = "6-NKT.rds")
dpi = 300
png(file="6-NKT-umap-split-sample.png", width = dpi*20, height = dpi*18, units = "px",res = dpi,type='cairo')
DimPlot(object = seu, reduction = 'umap',label = TRUE, split.by = 'TIG', ncol = 5)
dev.off()

png(file="6-NKT-umap-split-group.png", width = dpi*20, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = seu, reduction = 'umap',label = TRUE, split.by = 'tissue',ncol = 3)
dev.off()

png(file="6-NKT-umap-group-sample.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
DimPlot(object = seu, reduction = 'umap',label = TRUE, group.by = 'sample')
dev.off()










# DefaultAssay(cd8s) <- "RNA"
# DefaultAssay(seu) <- "RNA"
gs <- c("YBX3", "CD37", "TCL1A", "PLCG2", "KLF2", "CD79B", "CD74", "HLA-DRA", "HLA-DRB1", "LAPTM5", "BTG1", "HLA-DQA2", "HLA-DPB1", "CD69", "HLA-DQB1", "HLA-DMB", "LINC00926", "RCSD1", "STK17A", "HLA-DPA1")
# cd8s <- AddModuleScore(seu, list(gs), name = "Bnaive")
# VlnPlot(cd8s, "Bnaive1", pt.size = 0)
# gs <- c("LTB", "HLA-DPB1", "HLA-DRA", "HLA-DPA1", "MTRNR2L8", "GPR183", "HLA-DRB1", "CD52", "BANK1", "SMCHD1", "RPL27A", "MS4A1", "CD74", "HLA-A", "RPS20", "RPS12", "RPL13A", "RPL30", "ATP1B3", "TNFRSF13C")
# cd8s <- AddModuleScore(seu, list(gs), name = "Bmem_switched")
# VlnPlot(cd8s, "Bmem_switched1", pt.size = 0)
# gs <- c("BCL2A1", "DUSP2", "MYC", "MIR155HG", "CD83", "EGR3", "REL", "CD69", "NFKBID", "CCR7", "PLEK", "FABP5", "ZFP36L1", "HLA-DQA1", "NR4A1", "DDX21", "NFKB1", "HSP90AB1", "CD70", "PIKFYVE")
# cd8s <- AddModuleScore(seu, list(gs), name = "Bmem_unswitched")
# VlnPlot(cd8s, "Bmem_unswitched1", pt.size = 0)
# gs <- c("IGHA2", "IGHA1", "JCHAIN", "IGHM", "IGLC2", "IGLL5", "IGLC3", "SSR4", "IGKC", "DERL3", "ANKRD28", "VIM", "GADD45A", "MZB1", "ERN1", "IGHD", "HSPA5", "SEC11C", "NEAT1", "CYTOR")
# cd8s <- AddModuleScore(seu, list(gs), name = "Bplasma-lgA")
# VlnPlot(cd8s, "Bplasma-lgA1", pt.size = 0)
# gs <- c("IGHG1", "IGHG4", "IGHG3", "IGHG2", "IGHGP", "IGLC7", "IGKC", "MZB1", "PRDX4","DNAAF1", "IGLV3-1", "JSRP1", "FKBP11", "SSR4", "IGKV3-20", "XBP1", "ITM2C", "IGKV4-1", "PSAP", "NPC2")
# cd8s <- AddModuleScore(seu, list(gs), name = "Bplasma-lgG")
# VlnPlot(cd8s, "Bplasma-lgG1", pt.size = 0)+
# VlnPlot(cd8s, "Bplasma-lgA1", pt.size = 0)

#######02#######
#!/usr/bin/env Rscript

### title: DGE (Differential gene expression) in CD8+ T cells comparing expanded vs non-expanded TCRs
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
lib_loc="/home/data/refdir/Rlib"
library(devtools,lib.loc = lib_loc)

devtools::install_github("montilab/hypeR")
library(hypeR,lib.loc = lib_loc)
library(DropletUtils)
'%notin%' <- Negate('%in%')

colExp<-c("#006E82", "#AA0A3C")


### Compare all expanded CD8+ T cells vs non-expanded CD8+ T cells
type <- 'T_vs_N'

# Read-in integrated object
#seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, seu$integrated_snn_res.0.1 == '0')
seu <- ScaleData(seu)



#### DEG ####
Idents(seu) <- seu$tissue
markers <- FindMarkers(seu, ident.1 = 'Tumor', ident.2 = 'Normal', 
                       min.pct = 0, logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA')
markers$gene <- rownames(markers)

tops <- markers %>%
  filter(p_val_adj < 0.05) %>% 
  top_n(n = 10, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% 
                           filter(p_val_adj < 0.05) %>% 
                           top_n(n = -10, wt = avg_log2FC))
markers$diffexpressed <- ifelse(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, 'Up', 
                                ifelse(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, 
                                       'Down', 'NS'))
markers$label <- ifelse(markers$gene %in% tops$gene | markers$gene %in% tops$gene, 
                        markers$gene, NA)

ifelse(!dir.exists(file.path(paste0('~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B/', type))), 
       dir.create(file.path(paste0('~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B/', type)), recursive = T), FALSE)
write.csv(markers, paste0('~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B/', type, '/markers_MBM_sc_', type, '_clone_size.csv'), row.names = F)

# Plot
pdf(paste0('~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/B/', type, '/plots_volcano_MBM_sc_', type, '_clone_size.pdf'))
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + 
  geom_text_repel() + 
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) + 
  labs(color = 'Log2FC\n(Tumor vs\nNormal)') + 
  theme_bw() + 
  coord_cartesian(xlim = c(-3, 3)) + 
  ggtitle(paste0('Tumor vs Normal ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + 
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) + 
  labs(color = 'Log2FC\n(Tumor vs\nNormal)') + 
  theme_bw() + 
  coord_cartesian(xlim = c(-3, 3)) + 
  ggtitle(paste0('Tumor vs Normal ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
dev.off()


DefaultAssay(seu) <- 'RNA'   
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling',
                         'non-cycling')
# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, "markers.res.0.1.csv", row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>%
  ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap.pdf'), width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'ident', raster = T, assay = 'RNA')
dev.off()
# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(0, 4, 5),
                             'Activated B cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(7), 'Naïve B cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(1, 3, 6), 'Plasma cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(2), 'Doublets',
                             seu$cell_type_fine)
# Sa# Differential gene expression (DGE) based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_celltype.csv'), row.names = F)
# ve object and labels
saveRDS(seu, paste0(path.ct, 'data_', filename, '.rds'))
ct <- seu@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct, 'data_', filename, '_celltype.csv'))

# Differential gene expression (DGE) based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_celltype.csv'), row.names = F)
# Differential gene expression (DGE) based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_celltype.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>%
  ScaleData()
DefaultAssay(seu) <- 'integrated'
pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap_celltype.pdf'), width = 10,
    height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T,
          assay = 'RNA')
dev.off()
