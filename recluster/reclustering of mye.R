####前半部分是 ~/202241028/Lambrecht.TNBC/mye.Lambrechs.20241030.R  分群是Lambrechs
#####改了 ~/LAMbrechts_anti_PD1_chemo/0913tnbc/0913doublet/mye.Lambrechs.20241030.R 按照这里的RenameIdent
####后半部分是 T.mye??
library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
install.packages("NMF")
library(NMF)
library(registry)
library(rngtools)
library(purrr)
library(cowplot)
library(stringr)
# library(pkgmaker)
library(cluster)

seu <- subset(data.seurat.filtered, cell_type_fine %in% c('Myeloid cells','pDCs'))
CS10_YS <- seu
# CS10_YS[["percent.mt"]] <- PercentageFeatureSet(CS10_YS, pattern = "^MT-")
# CS10_YS <- subset(CS10_YS, subset = nFeature_RNA > 500 & percent.mt < 10)   
dim(CS10_YS)
seu.1 <- CS10_YS
mito.genes <- grep(pattern = "^MT-", x = rownames(seu.1), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(seu.1), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(seu.1), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(seu.1), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(seu.1), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(seu.1), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(seu.1), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(seu.1), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(seu.1), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(seu.1), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(seu.1), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(seu.1), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(seu.1), value = TRUE)
counts <- GetAssayData(seu.1, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
CS10_YS <- subset(seu.1, features = rownames(counts))
# dim(CS10_YS)
# > dim(seu)
# [1] 25288  5000
# > dim(CS10_YS)
# [1] 25146  5000


CS10_YS <- NormalizeData(object = CS10_YS)
CS10_YS <- FindVariableFeatures(object = CS10_YS)
dim(CS10_YS)

CS10_YS <- CellCycleScoring(
  object = CS10_YS,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

sf.list <- SplitObject(CS10_YS,split.by=c("patient_id"))
YS.anchors <- FindIntegrationAnchors(object.list = sf.list, dims = 1:30)
lapply(sf.list, dim)
YS.combined <- IntegrateData(anchorset = YS.anchors, dims = 1:30,k.weight=20)
DefaultAssay(YS.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
all.genes <- rownames(YS.combined)
YS.combined <- ScaleData(YS.combined, features = all.genes, vars.to.regress =c("nCount_RNA","percent.mt","S.Score","G2M.Score"))
YS.combined <- RunPCA(YS.combined, npcs = 50, features = VariableFeatures(object = YS.combined))

ElbowPlot(YS.combined)
dim(YS.combined@meta.data)
YS.combined <- RunUMAP(YS.combined, reduction = "pca", dims = 1:15)
YS.combined <- RunTSNE(YS.combined, reduction = "pca", dims = 1:15)
YS.combined <- FindNeighbors(YS.combined, reduction = "pca", dims = 1:15)

YS.combined <- FindClusters(YS.combined, resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))
YS.combined <- FindClusters(YS.combined, resolution = c(0.8))
table(YS.combined@active.ident)

clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(YS.combined@meta.data),value = TRUE),
                            function(x) length(unique(YS.combined@meta.data[,x])))
clusters_resolution
Idents(YS.combined) <- "integrated_snn_res.0.6"#!!!
DefaultAssay(YS.combined) <- "RNA"
YS.combined_markers <- FindAllMarkers( object = YS.combined,only.pos = T,
                                       test.use = "MAST")
YS.combined_markers <- YS.combined_markers[YS.combined_markers$p_val_adj < 0.01,]

harmonized_seurat2 <- RenameIdents(YS.combined,
                                     "11"="cDC1-CLEC9A",     # DC细胞群
                                     "3"="cDC2-CD1C",        # DC细胞群
                                     "12"="cDC3-LAMP3",      # DC细胞群
                                     "5"="pDC-LILRA4",      # DC细胞群
                                     "2"="Mono-FCN1",       # 单核细胞
                                     "4"="Macro-SPP1",      # 巨噬细胞
                                     "0"="Macro-C1QC",      # 巨噬细胞
                                     "1"="Macro-C1QC",      # 巨噬细胞
                                     "7"="Macro-FOLR2",     # 巨噬细胞
                                     "13"="Macro-FOLR2",    # 巨噬细胞
                                     "6"="Macro-CXCL10",    # 巨噬细胞
                                     "8"="T-B",          # 其他污染的细胞
                                     "10"="Fibroblast",     # 其他污染的细胞
                                     "9"="Epithelia"        # 其他污染的细胞
  )


harmonized_seurat2$cell_type_fine <- Idents(harmonized_seurat2)
DimPlot(harmonized_seurat2,label = T)
unique(harmonized_seurat2$cell_type_fine)
harmonized_seurat2$barcode <- colnames(harmonized_seurat2)
celltype <- harmonized_seurat2@meta.data %>% select('barcode', 'cell_type_fine')
write.csv(celltype, 'mye_celltype.csv')

Mye <- readRDS("/home/data/t180330/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Myeloid cells/int.Mye.RENAMEDIDENT.RDS")
unique(Mye$seurat_clusters)
levels(Mye$seurat_clusters)
# 假设seu对象是您的Seurat对象，且它的细胞分类信息存储在seurat_clusters这个metadata中
# 确定要排除的细胞类别的名称
'%notin%' <- Negate('%in%')
excluded_clusters <- c("T-B", "Fibroblast", "Epithelia")
# 使用subset函数排除这些细胞，这里直接使用!符号表示不在这些类别中
Mye.filt <- subset(Mye, seurat_clusters %notin% excluded_clusters)
dim(Mye.filt)
# > dim(Mye.filt)
# [1] 30959  3497


Mye.filt <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Myeloid cells/co-inhi-co-stimu/int.Mye.RENAMEDIDENT.filt.RDS")
unique(Idents(Mye.filt))
DimPlot(Mye.filt,label = T)
saveRDS(Mye.filt,"Mye.filt.RDS")
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
logfc.threshold = 0.25,test.use = "MAST") #
write.csv(markers, 'markers_lungs_all_cells.csv', row.names = F)
DefaultAssay(Mye.filt) <- 'RNA'
Mye.filt <- CellCycleScoring(Mye.filt, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(Mye.filt) <- 'integrated'
Mye.filt$cell_cycle <- ifelse(Mye.filt$G2M.Score > 0.13 | Mye.filt$S.Score > 0.13, 'cycling', 'non-cycling')

DefaultAssay(Mye.filt) <- 'RNA'
Mye.filt <- CellCycleScoring(Mye.filt, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(Mye.filt) <- 'integrated'
Mye.filt$cell_cycle <- ifelse(Mye.filt$G2M.Score > 0.13 | Mye.filt$S.Score > 0.13, 'cycling', 'non-cycling')

unique(Mye.filt$seurat_clusters)

# Fine cell type assignment after manual annotation based on DGE
Mye.filt$cell_type_fine <- ifelse(Mye.filt$seurat_clusters %in% c("Macro-C1QC","Macro-FOLR2","Macro-CXCL10","Macro-SPP1"), 'Macrophages', 'NA')
Mye.filt$cell_type_fine <- ifelse(Mye.filt$seurat_clusters %in% c("Mono-FCN1"), 'Monocytes', Mye.filt$cell_type_fine)
# Mye.filt$cell_type_fine <- ifelse(Mye.filt$seurat_clusters %in% c(3, 16), 'remove', Mye.filt$cell_type_fine)
Mye.filt$cell_type_fine <- ifelse(Mye.filt$seurat_clusters %in% c("cDC1-CLEC9A", "cDC3-LAMP3", "cDC2-CD1C","pDC-LILRA4"), 'DCs', Mye.filt$cell_type_fine)
Mye.filt$barcode <- colnames(Mye.filt)
saveRDS(Mye.filt, 'Mye.filt.cell_type_fine.rds')

celltype <- Mye.filt@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, "myeloid_celltype.csv")

Idents(harmonized_seurat2) <- harmonized_seurat2$seurat_clusters
YS.combined_markers <- FindAllMarkers(object = harmonized_seurat2,only.pos =TRUE,
                                                                               test.use = "MAST")
head(YS.combined_markers)
library(dplyr)
# cluster.averages_Pan02_CD8 <- AverageExpression(harmonized_seurat2, return.seurat = TRUE)
top30 <- YS.combined_markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
cell_colors <- c(
  "#BC6955",
  "#F0C6BE",
  "#F9E7C0",
  "#DFD57C",
  "#ACA7C2",
  "#CAE5D7",
  "#A6C694",
  "#CDE0E8",
  "#42A6E5",
  "#AF5367",
  "#61C1C5",
  "#ACA442"
)
harmonized_seurat2 <- ScaleData(harmonized_seurat2, features = rownames(harmonized_seurat2))
# ph=DoHeatmap(harmonized_seurat2, features = top30$gene, group.colors=cell_colors ,hjust=F, angle = 45, size = 4) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))   +
#   theme (axis.text.y = element_blank()) + NoLegend()
# ph
ph=DoHeatmap(harmonized_seurat2, features = top30$gene, group.colors=cell_colors ,hjust=F,draw.lines = T, angle = 45, size = 4) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))   +
  theme (axis.text.y = element_blank()) 
ph
dev.off()
write.table(top30,"alltop30.hm.marker.tsv",sep = '\t')
write.table(YS.combined_markers,"all.mye.hm.marker.tsv",sep = '\t')
