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

Tcell <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/T.rds")
DimPlot(Tcell,label = T)
DimPlot(Tcell,label = T,repel=F)
unique(Idents(Tcell))
# 定义细胞身份列表
cell_identities <- c(
  "CD4/CD8-PLCG2", "CD8-CXCL13", "CD8-GZMK", "CD8-CCL5", "CD4-CCR7",
  "NK-FGFBP2", "CD4-CXCL13", "Epithelia", "CD8-gdT-TRDC", "CD4-FOXP3",
  "Cycling cells", "NK-FCGR3A", "Stress cells"
)
# 要移除的身份
idents_to_remove <- c("CD4/CD8-PLCG2", "Epithelia")
# 使用 setdiff 函数移除指定的身份
substituted_identities <- setdiff(cell_identities, idents_to_remove)
# 打印更新后的列表
print(substituted_identities)
Tcell$seurat_clusters <- Idents(Tcell)

T.filt <- subset(Tcell,seurat_clusters %in% substituted_identities )
DimPlot(T.filt,label = T)
saveRDS(T.filt,"T.filt.RDS")



# Cell cycle assignment
DefaultAssay(T.filt) <- 'RNA'
T.filt <- CellCycleScoring(T.filt, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(T.filt) <- 'integrated'
T.filt$cell_cycle <- ifelse(T.filt$G2M.Score > 0.15 | T.filt$S.Score > 0.1, 'cycling', 'non-cycling')
unique(T.filt$seurat_clusters)
T.filt$cell_type_fine <- T.filt$seurat_clusters




T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("CD4-FOXP3"), 'Tregs', 'NA')
T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("CD4-CCR7", "CD4-CXCL13"), 'CD4+ T cells', T.filt$cell_type_fine)
T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("CD8-CXCL13","CD8-GZMK","CD8-CCL5","CD8-gdT-TRDC"), 'CD8+ T cells', T.filt$cell_type_fine)
T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("NK-FCGR3A", "NK-FGFBP2"), 'NK cells', T.filt$cell_type_fine)
T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("Cycling cells"), 'Cycling NK/T cells', T.filt$cell_type_fine)
T.filt$cell_type_fine <- ifelse(T.filt$seurat_clusters %in% c("Stress cells"), 'Stressed NK/T cells', T.filt$cell_type_fine)

DimPlot(T.filt,group.by = "cell_type_fine")

# Save object and labels
T.filt$barcode <- colnames(T.filt)
saveRDS(T.filt, 'T.filt.cell_type_fine.rds')
library("dplyr")
celltype <- T.filt@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'tcells_celltype.csv')
