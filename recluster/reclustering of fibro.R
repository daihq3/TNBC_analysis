# ~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/DC/fib.NT.R
Mye <- subset(CombinedMFLN920, cell_type_main %in% c("Fibroblasts","Pericytes","Endothelium"))
# [1] 25288   354
dim(Mye)
# rm(DCs)
table(CombinedMFLN920$cell_type_main)
# rm(Combined_MFLN9_20)
# [1] 25288  5000
# data.seurat.filtered <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/seu.withcluster-rename.rds")
# data.seurat.filtered <- Mye
saveRDS(Mye,"data.seurat.filtered.rds")
setwd("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts")
# Mye <- readRDS("~/202241028/Lambrecht.TNBC/data.seurat.filtered.rds")
# data.seurat.filtered <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/seu.withcluster-rename.rds")
seu <- Mye
# rm(data.seurat.filtered)
DefaultAssay(seu) <- "RNA"

CS10_YS <- seu
CS10_YS[["percent.mt"]] <- PercentageFeatureSet(CS10_YS, pattern = "^MT-")
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
# > dim(CS10_YS)
# [1] 30959  3823

CS10_YS <- CellCycleScoring(
  object = CS10_YS,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
dim(YS.combined.filt)

# CS10_YS$
sf.list <- SplitObject(CS10_YS,split.by="TIG")
YS.anchors <- FindIntegrationAnchors(object.list = sf.list, dims = 1:30)
lapply(sf.list, dim)
YS.combined <- IntegrateData(anchorset = YS.anchors, dims = 1:30)
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

YS.combined <- FindClusters(YS.combined, resolution = c(0.3,0.4,0.5,0.6,0.7))
YS.combined <- FindClusters(YS.combined, resolution = c(0.1,0.2))
table(YS.combined@active.ident)
dim()
clusters_resolution<-sapply(grep("^integrated_snn_res",colnames(YS.combined@meta.data),value = TRUE),
                            function(x) length(unique(YS.combined@meta.data[,x])))
clusters_resolution
Idents(YS.combined) <- "integrated_snn_res.0.2"#!!!
DefaultAssay(YS.combined) <- "RNA"
Idents(YS.combined)
YS.combined_markers <- FindAllMarkers( object = YS.combined,only.pos = T,
                                       test.use = "MAST")
YS.combined_markers2 <- FindAllMarkers( object = harmonized_seurat2,only.pos = T,
                                       test.use = "wilcox")
YS.combined_markers <- YS.combined_markers[YS.combined_markers$p_val_adj < 0.01,]
write.table(YS.combined_markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
# YS.combined_markers <- YS.combined_markers %>% group_by(gene) %>% top_n(n = 1, wt = avg_logFC)
write.csv(YS.combined_markers,"~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts/Res.0.3.9cluster.of.stromal.cells")
write.csv(YS.combined_markers,"~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts/Res.Stromal.7cluster.of.stromal.cells.csv")
write.csv(YS.combined_markers,"~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts/Res.Stromal.7cluster.of.stromal.cells.csv")
write.csv( ,"~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts/top_stroma_markers.csv")
# 然后，将ident列转换为因子
Ery.integrated.markers <- FindAllMarkers( object = harmonized_seurat2,only.pos = T,logfc.threshold = 0.25,
                                          test.use = "wilcox")



top_markers = YS.combined_markers %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  pull(gene) %>% as.character()
Idents(YS.combined) <- "integrated_snn_res.0.2"#!!!

YS.combined <- RenameIdents(YS.combined, 
                                                                     `0` = "Fibro-PI16",
                                                                     `1` = "Fibro-CXCL1",
                                                                     `2` = "Fibro-POSTN",
                                                                     `3` = "Endo-ACKR1",
                                                                     `4` = "Peri-RGS5",
                                                                     `5` = "Endo-RGCC",
                                                                     `6` = "Doublet",
                                                                     `7` = "Doublet",
                                                                     `8` = "Endo-LYVE1")
YS.combined <- RenameIdents(YS.combined, 
                                                                     `0` = "Fibro-PI16",
                                                                     `1` = "Fibro-CXCL1",
                                                                     `2` = "Fibro-POSTN",
                                                                     `3` = "Endo-ACKR1",
                                                                     `4` = "Peri-RGS5",
                                                                     `5` = "Endo-RGCC",
                                                                     `6` = "Doublet",
                                                                     `7` = "Doublet",
                                                                     `8` = "Endo-LYVE1")
# 加载dplyr包
library(dplyr)
Idents(YS.combined)
# 假设YS.combined_markers是FindAllMarkers函数返回的标记表格
# 使用filter函数去除cluster为"Doublet"的行
YS.combined_markers_clean <- YS.combined_markers %>%
  filter(cluster != "Doublet")

head(YS.combined_markers)
library(dplyr)
top30 <- YS.combined_markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)

# 然后，将ident列转换为因子
Idents(YS.combined.filt) <- factor(Idents(YS.combined.filt), 
                                           levels = c("Fibro-PI16", "Fibro-CXCL1", "Fibro-POSTN", 
                                                      "Endo-ACKR1","Endo-RGCC", "Endo-LYVE1", "Peri-RGS5"))


DimPlot(YS.combined.filt,cols = col_fibro)
                            
top_basal_markers <- YS.combined_markers_clean %>%
  group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP|^AC|^AL|^LINC")) %>%
  dplyr::slice_max(n = 30, order_by = avg_log2FC) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  as.character()
p1 = DoHeatmap(YS.combined.filt, features = top_basal_markers, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2,group.colors = col_fibro,draw.lines = F, raster = F)  +
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_fibro_allcells_v3.pdf", width = 7, height = 8)

# Load the dplyr package if you haven't already
library(dplyr)

# Count the number of genes per cluster
YS.combined_markersclean2 %>%
  group_by(cluster) %>%
  summarise(number_of_genes = n())

DefaultAssay(YS.combined.filt) <- "RNA"
YS.combined.filt <- ScaleData(YS.combined.filt, features = rownames(YS.combined.filt))
ph=DoHeatmap(YS.combined.filt, features = top_basal_markers, group.colors= col_fibro ,hjust=F,draw.lines = F, angle = 45, size = 4) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 4, name = "RdYlBu"))) + theme(axis.text.y = element_text(size=2.5), axis.text.x = element_text(size=6.5))

theme (axis.text.y = element_blank()) 
ph
ggsave(plot = ph, filename = "heatmap_cellstate_fibro_allcells_v4.pdf", width = 7, height = 8)
# Access the data from the plot object
plot_data <- ggplot_build(p1)$data[[1]]

# Extract gene names from the plot data
gene_names_from_plot <- plot_data$label
print(gene_names_from_plot)

dim(YS.combined)
# [1] 30959 13773
'%notin%' <- Negate('%in%')
YS.combined.filt <- subset(YS.combined, integrated_snn_res.0.2 %notin% c(6,7))


dim(YS.combined.filt)
# [1] 30959 12625
# YS.combined <- ScaleData(YS.combined)
# DoHeatmap(YS.combined, top_markers)#, group.colors = col_fibro)
fibro <- YS.combined.filt
DefaultAssay(fibro) <- "integrated"
YS.combined_markersclean2 <- YS.combined_markers2 %>%
  filter(cluster != "Doublet")

tum_sub1 = AverageExpression(fibro, features = unique(YS.combined_markersclean2$gene), return.seurat = T, assays = "RNA")
# tum_sub1 = AverageExpression(fibro, features = unique(top_basal_markers), return.seurat = T, assays = "RNA")
Idents(tum_sub1)
DefaultAssay(tum_sub1) <- "RNA"
tum_sub1 <- ScaleData(tum_sub1)

DefaultAssay(tum_sub1) <- "RNA"

tum_sub1 <- ScaleData(tum_sub1, features = rownames(tum_sub1))

ph=DoHeatmap(tum_sub1, features = YS.combined_markersclean2$gene, group.colors=col_fibro ,hjust=F,draw.lines = F, angle = 45, size = 4) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))   +
  theme (axis.text.y = element_blank()) 
ph

# 代跑1 ---------------------------------------------------------------------
col_fibro = c("#A5D1B0", "#CE8ABD", "#FFF7C1", "#ADD3F4", "#F7C9CF", "#7CA3BB", "#5678D6")
consistentcolors2 <- c("#006E82", "#AA0A3C")
p1 = DimPlot(harmonized_seurat2,cols = col_fibro)
ggsave(plot = p1, filename = "DimPlot.stromal.pdf", width = 6, height = 4)
p2 = DimPlot(harmonized_seurat2,cols = col_fibro,split.by = "tissue")
p4 = DimPlot(harmonized_seurat2,cols = consistentcolors2,group.by = "tissue")
# ggsave(plot = p2, filename = "DimPlot.stromal.pdf", width = 10, height = 4)
ggsave(plot = p2, filename = "DimPlot.stromal.split.byTNC.pdf", width = 10, height = 4)
ggsave(plot = p4, filename = "DimPlot.stromal,group.byTNC.pdf", width = 5, height = 4)

# 代跑2 ---------------------------------------------------------------------

p1 = DoHeatmap(tum_sub1, features = unique(YS.combined_markersclean2$gene) , assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2,group.colors = col_fibro,draw.lines = F, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_fibro_allcells_v3.pdf", width = 4, height = 6)
p1 = DoHeatmap(tum_sub1, features = unique(YS.combined_markersclean2$gene) , assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2,group.colors = col_fibro,draw.lines = F, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=1.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_fibro_allcells_v3.1.pdf", width = 4, height = 6)
col_fibro = c("#A5D1B0", "#CE8ABD", "#FFF7C1", "#ADD3F4", "#F7C9CF", "#7CA3BB", "#5678D6")
p1 = DoHeatmap(harmonized_seurat2, features = top_basal_markers, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2,group.colors = col_fibro,draw.lines = F, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_fibro_allcells_v3.pdf", width = 4, height = 10)
fibro = AddModuleScore(fibro,features = list(vec_chemokines), assay = "RNA", name = "chemokines")
v1=VlnPlot(fibro, features = "chemokines1", pt.size = 0, group.by = "final_group", cols = col_fibro)
ggsave(plot = v1, filename = "vlnplot_chemokines_v1.pdf", width = 4, height = 3)

FeaturePlot(sc_hc,c("Aoc3","Acta2"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
FeaturePlot(sc_hc,c("Sox6","Wnt5a"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
FeaturePlot(sc_hc,c("Maff","Cd55"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
FeaturePlot(sc_hc,c("Ier3","Ptgs2"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
FeaturePlot(sc_hc,c("Cxcl12","Adamdec1"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
FeaturePlot(sc_hc,c("Edil3","Hapln1"), cols.use = c("lightgrey","blue"), pt.size = 0.6)
marker_list <- c("Epcam","Krt19","Kit","Ano1","Myh11","Des",
                 "Rgs5","Pdgfrb","Pecam1","Lyve1",
                 "S100b","Gfap","Cd52","Ptprc",
                 "Acta2","Dpt","Col6a2","Col1a2")
marker_list <- c("Epcam","Krt19","Kit","Ano1","Myh11","Des",
                 "Rgs5","Pdgfrb","Pecam1","Lyve1",
                 "S100b","Gfap","Cd52","Ptprc",
                 "Acta2","Dpt","Col6a2","Col1a2",)

gs=c("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4,PROX1")
gs <- strsplit("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4,PROX1,PLVAP,MKI67,APOE,FRZB,FBLN1", split = ",")[[1]]
print(gs)
p <- DotPlot(object = YS.combined.filt,features = marker_list_upper, assay = "RNA", scale = TRUE) +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20))
p
gs <- strsplit("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4,PROX1,PLVAP,MKI67,APOE,FRZB,FBLN1", split = ",")[[1]]



gs <- strsplit("PI16,CXCL1,POSTN,ACKR1,RGS5,RGCC,LYVE1", split = ",")[[1]]
gs <- strsplit("PI16,CXCL1,POSTN,ACKR1,RGCC,LYVE1,RGS5", split = ",")[[1]]

harmonized_seurat2 <- YS.combined.filt
harmonized_seurat2@meta.data$cell_type_fine <- fct_relevel(harmonized_seurat2@meta.data$cell_type_fine, major_celltype_order)
p <- DotPlot(object = harmonized_seurat2,features = gs, assay = "RNA", scale = TRUE,group.by = "cell_type_fine") +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12))
# pdf(file="tcr-heatmap-seurat.pdf", width = 9, height = 3.5)
# pp = DotPlot(YS.combined.filt, features = gs,cols = c('white', "#308a46"),dot.scale=4.5,group.by = "cell_type_fine") + RotatedAxis()
# pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
#   guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
#   theme(axis.line = element_line(size = 1))
# print(pp)
# dev.off()
YS.combined.filt
YS.combined.filt$cell_type_fine <- Idents(YS.combined.filt)
major_celltype_order <- rev(levels(YS.combined.filt$cell_type_fine))
major_celltype_order
library(dplyr)
library(forcats)


nCoV_markers_rev <- gs
pdf(file="tcr-heatmap-seurat.pdf", width = 9, height = 3.5)
pp = DotPlot(harmonized_seurat2, features = nCoV_markers_rev,cols = c('white', "#308a46"),dot.scale=4.5,group.by = "cell_type_fine") + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 1))
print(pp)
dev.off()

gs

write.csv(YS.combined_markers2,"YS.combined_markers2.csv")
# 加载dplyr包
library(dplyr)

# 假设YS.combined_markers2是你的数据框，且包含'gene'和'cluster'两列

# 计算每个簇特有的基因
unique_genes_per_cluster <- YS.combined_markers2 %>%
  group_by(cluster) %>%
  # 获取当前簇的基因
  summarise(genes_in_cluster = list(gene)) %>%
  # 检查每个基因是否只出现在当前簇
  mutate(unique_genes = lapply(genes_in_cluster, function(genes) {
    # 获取所有簇的基因
    all_genes <- YS.combined_markers2$gene
    # 找出只出现在当前簇的基因
    genes[!(genes %in% all_genes[YS.combined_markers2$cluster != cluster])]
  })) %>%
  # 提取簇名称和特有基因列表
  select(cluster, unique_genes)

# 查看结果
print(unique_genes_per_cluster)

# For this annotation, we used known markers, including GJA4 for arterial BECs, CA4 for capillary BECs and ACKR1 for venous BECs29,30,31 (Fig. 2c,d and Supplementary Table 4).
# we used known markers, including GJA4 for arterial BECs,

gs=c("PECAM1,JAM2,PROX1,DCN,ACTA2,GJA4,CA4,ACKR1")
gs=c("PECAM1,JAM2,PROX1,DCN,ACTA2,GJA4,CA4,ACKR1,RGCC")
gs=c("PDGFRA,COL1A1,DCN,ACTA2,PDGFRB,MCAM,PECAM1,CD34,VWF")


gs <- strsplit("PDGFRA,COL1A1,DCN,PECAM1,CD34,VWF,ACTA2,PDGFRB,MCAM", split = ",")[[1]]

     #,JAM2,ACTA2,DCN,ACTA2,GJA4,CA4,ACKR1,RGCC")


gs=c("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4")
gs=strsplit(gs)[[1]]
gs=strsplit("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4")[[1]]
gs <- strsplit("PECAM1,VWF,SOX17,GJA4,ACKR1,SELP,RGCC,CA4", split = ",")[[1]]

FeaturePlot(tum_sub1,features = c())

marker_list_upper <- toupper(marker_list)
print(marker_list_upper)

marker_spec <- c("Epithelial","Epithelial","ICC","ICC","Smooth muscle",
                 "Smooth muscle","Pericyte","Pericyte","Endothelial","Lymphatic",
                 "Glial","Glial","Immune","Immune","Fibroblast",
                 "Fibroblast","Fibroblast","Fibroblast")

DimPlot(YS.combined,label = T)


####~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/Fibroblasts/stromal.hm.ANGIO.IMME.R
Idents(YS.combined) <- factor(Idents(YS.combined),
                                   levels = c("Fibro-PI16", "Fibro-CXCL1", "Fibro-POSTN",
                                              "Endo-ACKR1","Endo-RGCC", "Endo-LYVE1", "Peri-RGS5","Doublet"))
YS.combined_markers <- FindAllMarkers( object = YS.combined,only.pos = T,
                                       test.use = "MAST")
YS.combined_markersclean2 <- YS.combined_markers %>%
  filter(cluster != "Doublet")
YS.combined_markers33 <- YS.combined_markersclean2[YS.combined_markersclean2$avg_log2FC > 2,]
YS.combined_markers2 <- YS.combined_markers33[YS.combined_markers33$p_val_adj < 0.01,]
write.csv(YS.combined_markers2,"YS.combined_markers2.csv")
fibro <- YS.combined.filt
DefaultAssay(fibro) <- "integrated"
tum_sub1 = AverageExpression(fibro, features = unique(YS.combined_markersclean2$gene), return.seurat = T, assays = "RNA")
DefaultAssay(tum_sub1) <- "RNA"
tum_sub1 <- ScaleData(tum_sub1)
View(top30)
col_fibro = c("#A5D1B0", "#CE8ABD", "#FFF7C1", "#ADD3F4", "#F7C9CF", "#7CA3BB", "#5678D6")
ph=DoHeatmap(tum_sub1, features = YS_combined_markers2_1$gene, group.colors=col_fibro ,hjust=F,draw.lines = F, angle = 45, size = 4) +   scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") +
    theme (axis.text.y = element_blank())
ph
ggsave(plot = ph, filename = "heatmap_cellstate_fibro_allcells_v3.pdf", width = 4, height = 6)

gs <- strsplit("PDGFRA,COL1A1,DCN,PECAM1,CD34,VWF,ACTA2,PDGFRB,MCAM", split = ",")[[1]]

p <- DotPlot(object = harmonized_seurat2,features = gs, assay = "RNA", scale = TRUE,group.by = "cell_type_fine") +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12))
ggsave(plot = p, filename = "Pi16_Dotplot_cluster.pdf", width = 6.4, height = 4)

# genes<-c("APOE", "S100A8", "FTL3", "LILRA4")
gs <- strsplit("PDGFRA,COL1A1,PECAM1,VWF,,ACTA2", split = ",")[[1]]
p3 <- FeaturePlot(harmonized_seurat2, features = gs,order=T,
                  reduction = "umap",pt.size = 1,  cols = c("grey", "red"),combine = FALSE)
p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + NoLegend()))
CombinePlots(plots = p3,ncol=2)
ggsave("gene1s.pdf",height =10,width =9)


# 定义正确的顺序
correct_order <- c("Fibro-PI16", "Fibro-CXCL1", "Fibro-POSTN", "Endo-ACKR1", "Endo-RGCC", "Endo-LYVE1", "Peri-RGS5")

CellsClusters <- data.frame(Cell = names(Idents(Endo)), 
                            CellType =  (Idents(Endo)), 
                            stringsAsFactors = FALSE) 
head(CellsClusters)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pbmc <- Endo
pbmc <- progeny(pbmc, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)  #"Human" Mouse
pbmc@assays$progeny 
pbmc@assays$progeny %>%dim()
pbmc@assays$progeny@data[,1:19]
library(viridis)
library(dplyr)
library(tidyr)
library(tibble)
## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
pbmc <- Seurat::ScaleData(pbmc, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(pbmc, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 14, treeheight_row = 0,
                        color=myColor, breaks = progenyBreaks, 
                        main = "", angle_col = 90,
                        treeheight_col = 0,  border_color = "black")
# 创建一个包含正确顺序的向量
order_vector <- c("Fibro-PI16", "Fibro-CXCL1", "Fibro-POSTN", "Endo-ACKR1", "Endo-RGCC", "Endo-LYVE1", "Peri-RGS5")

# 使用arrange函数根据顺序向量来排列数据框
YS.combined_markers2 <- YS.combined_markers2 %>%
  mutate(cluster_order = match(cluster, order_vector)) %>%
  arrange(cluster_order) %>%
  select(-cluster_order)

# 查看结果
head(YS.combined_markers2)


Mye.filt <- harmonized_seurat2
meta_data <- Mye.filt@meta.data  # 提取元数据
cell_state <- Idents(Mye.filt)  # 提取细胞状态（聚类）
meta_data$Cell_state <- cell_state  # 将细胞状态添加到元数据中
# 假设meta_data中包含器官和组织信息
# 请检查meta_data是否包含 `Organ_orig` 和 `Tissue` 两个列，如果没有，请替换为正确的列名
crosstab <- Mye.filt@meta.data
summary<-table(meta_data[,c('Cell_state','tissue')])#tissuegroup在数据均一致，不同于是否缩写
roe<-as.data.frame(ROIE(summary))
library("pheatmap")
pheatmap(roe,display_numbers=TRUE,number_color="black",cluster_row=FALSE,cluster_col=FALSE,
         color=colorRampPalette(c("#E4F0F0","#F33D0C"))(50))#画图

VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
Endo <- subset(harmonized_seurat2,cell_type_fine %in% c("Endo-ACKR1", "Endo-RGCC", "Endo-LYVE1"))
VlnPlot(Endo,features = c("ACTA2"), pt.size = 0, ncol = 1,cols = fibcol) + NoLegend()

VlnPlot(Endo,features = c("KCNE3","ANGPT2"))
VlnPlot(Endo,features = c("GJA4"))
VlnPlot(Endo,features = c("RAMP3"))
VlnPlot(Endo,features = c("ACTA2"),cols =  c("#ADD3F4", "#F7C9CF", "#7CA3BB"), pt.size = 0, ncol = 1) + NoLegend()
Endo <- subset(harmonized_seurat2,cell_type_fine %in% c("Fibro-PI16", "Fibro-CXCL1", "Fibro-POSTN"))
VlnPlot(Endo,features = c("CXCL1","CXCL3","CXCL6","LIF"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo,features = c("LRRC15"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo,features = c("GJA4"))
VlnPlot(Endo,features = c("RAMP3"))
VlnPlot(Endo,features = c("ACTA2"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 1) + NoLegend()
VlnPlot(Endo,features = c("PI16","CD34","DPT"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("COL1A2", "POSTN","SPP1"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("COL1A2", "POSTN","SPP1"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("ADM", "ANKRD37", "C1QTNF6", "CLDN5", "COL4A1", "COL4A2", "COTL1", "DLL4", "EDNRB"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", "INHBB", "JUP", "KCNE3", "KCNJ8", "LAMA4"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c( "LAMB1", "LXN", "MARCKSL1", "MCAM", "MEST", "N4BP3", "NID2", "NOTCH4", "PLOD1"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("PLXND1", "PMEPA1", "PTN", "RAMP3", "RBP1", "RGCC", "RHOC", "TRP53ILL"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
VlnPlot(Endo,features = c("CD36"),cols =  c("#A5D1B0", "#CE8ABD", "#FFF7C1"), pt.size = 0, ncol = 3) + NoLegend()
fibcol=c("#A5D1B0", "#CE8ABD", "#FFF7C1")
endcol=c("#ADD3F4", "#F7C9CF", "#7CA3BB")
# progeny
# used PROGENy to estimate the activity

result <- YS.combined_markers_clean %>%
  filter(cluster == "Fibro-POSTN") %>%
  arrange(desc(avg_log2FC)) %>%
  select(gene, avg_log2FC)

# 查看结果
print(result)

write.table(result,file="expression.rnk",quote=F,sep="\t",row.names=F,co)
YS.combined_markers <- FindAllMarkers( object = YS.combined,
                                       test.use = "MAST")

result <- all.Tcell.hm.marker %>%
  filter(cluster == "CD8-CXCL13") %>%
  arrange(desc(avg_log2FC)) %>%
  select(gene, avg_log2FC)

# 查看结果
print(result)

write.table(result,file="expression.rnk",quote=F,sep="\t",row.names=F,co)
VlnPlot(Endo, features = c("HSPG2", "ANGPT2", "HIF1A", "MMP2", "CTGF", "NOTCH1"))
VlnPlot(Endo, features = c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1"))
VlnPlot(Endo, features = c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1"))
VlnPlot(Endo, features = c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1"))
VlnPlot(Endo, features = c("COL1A2, POSTN, FAP, and ACTA2"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("COL1A1", "POSTN","FAP", "ACTA2"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("CXCL1", "IL6","CCL2", "IL1R1"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("CCND1", "COL8A1","MYLK", "CD73"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("LIF", "IL6","CXCL1", "CXCL3"), pt.size = 0, ncol = 1,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("MYH11", "ACTA2", "MCAM", "CAV1"), pt.size = 0, ncol = 2,cols = fibcol) + NoLegend()
VlnPlot(Endo, features = c("CTHRC1", "ACTA2", "POSTN", "INHBA","SFRP2","SULF1","COL8A1"),   = 0, ncol = 2,cols = fibcol) + NoLegend()

# # 
# # including CTHRC1, POSTN, INHBA, SFRP2,
# # SULF1, COL8A1, COL11A1 and LRRC15)
# library(Seurat)
# library(ggplot2)
# modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
#   p <- VlnPlot(obj, features = feature, pt.size = pt.size, ...) +
#     xlab("") +
#     ylab(feature) +
#     ggtitle("") +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_text(), # 确保 y 轴文本是显示的
#       axis.ticks.x = element_blank(),
#       axis.ticks.y = element_line(),
#       axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
#       plot.margin = plot.margin
#     ) +
#     coord_cartesian(ylim = c(0, 6)) # 设置 y 轴刻度范围为 0 到 6
#   
#   return(p)
# }
# 
# ## main function
# StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
#   plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
#   plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
#     theme(axis.text.x=element_text(), axis.ticks.x = element_line())
#   p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
#   return(p)
# }
# 
# 
# 
# StackedVlnPlot(Endo, c("HSPG2", "ANGPT2", "HIF1A", "MMP2", "CTGF", "NOTCH1"), pt.size=0, cols=my36colors)
# StackedVlnPlot(Endo, c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1"), pt.size=0, cols=my36colors)

#配色方案
my36colors <- endcol

# Endo$cell_type_fine
vln.dat=FetchData(Endo,c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1","cell_type_fine","tissue"))
consistentcolors2 <- c("#006E82", "#AA0A3C")
vln.dat$Cell <- rownames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","cell_type_fine","tissue"), 
                               measure.vars = c("CCL2", "BIRC3","HES1", "HLA-DQB1", "HLA-DQA1", "ICAM1"),
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(cell_type_fine,tissue, gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值



ggplot(vln.dat.melt, aes(x = cell_type_fine,y = Expr, fill = tissue))+
  
  geom_split_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_fill_manual(values = c("#006E82", "#AA0A3C")) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 12,color = "black"), axis.text.y=element_text(size = 12,color = "black"),axis.ticks.x = element_line(),panel.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  
  stat_compare_means(data = vln.dat.melt, aes(x = cell_type_fine,y = Expr),
                     
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      
                                      symbols = c("***", "**", "*", "-")),
                     method = "wilcox",
                     
                     label = "p.signif",
                     
                     label.y = max(vln.dat.melt$Values),
                     
                     hide.ns = F)+
  
  theme_classic()+
  
  xlab("")+
  
  ylab("")+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        
        legend.position = "top",
        
        #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
        
        legend.justification = "right")



vln.dat=FetchData(Endo,c("HSPG2", "ANGPT2", "HIF1A", "MMP2", "CTGF", "NOTCH1","cell_type_fine","tissue"))
consistentcolors2 <- c("#006E82", "#AA0A3C")
vln.dat$Cell <- rownames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","cell_type_fine","tissue"), 
                               measure.vars = c("HSPG2", "ANGPT2", "HIF1A", "MMP2", "CTGF", "NOTCH1"),
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(cell_type_fine,tissue, gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值



ggplot(vln.dat.melt, aes(x = cell_type_fine,y = Expr, fill = tissue))+
  
  geom_split_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_fill_manual(values = c("#006E82", "#AA0A3C")) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 12,color = "black"), axis.text.y=element_text(size = 12,color = "black"),axis.ticks.x = element_line(),panel.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  
  stat_compare_means(data = vln.dat.melt, aes(x = cell_type_fine,y = Expr),
                     
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      
                                      symbols = c("***", "**", "*", "-")),
                     method = "wilcox",
                     label = "p.signif",
                     
                     label.y = max(vln.dat.melt$Values),
                     
                     hide.ns = F)+
  
  theme_classic()+
  
  xlab("")+
  
  ylab("")+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        
        legend.position = "top",
        
        #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
        
        legend.justification = "right")
vln.dat=FetchData(Endo,c("HSPG2", "SERPINE1", "COL4A2", "LAMA4", "COL4A1", "SPARC","cell_type_fine","tissue"))
consistentcolors2 <- c("#006E82", "#AA0A3C")
vln.dat$Cell <- rownames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","cell_type_fine","tissue"), 
                               measure.vars = c("HSPG2", "SERPINE1", "COL4A2", "LAMA4", "COL4A1", "SPARC"),
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(cell_type_fine,tissue, gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值



ggplot(vln.dat.melt, aes(x = cell_type_fine,y = Expr, fill = tissue))+
  
  geom_split_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_fill_manual(values = c("#006E82", "#AA0A3C")) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 12,color = "black"), axis.text.y=element_text(size = 12,color = "black"),axis.ticks.x = element_line(),panel.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  
  stat_compare_means(data = vln.dat.melt, aes(x = cell_type_fine,y = Expr),
                     
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      
                                      symbols = c("***", "**", "*", "-")),
                     method = "wilcox",
                     label = "p.signif",
                     
                     label.y = max(vln.dat.melt$Values),
                     
                     hide.ns = F)+
  
  theme_classic()+
  
  xlab("")+
  
  ylab("")+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        
        legend.position = "top",
        
        #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
        
        legend.justification = "right")

# sce <- as.SingleCellExperiment(Endo)
# PTtrajectory_A <- Endo
# # Add metadata from Seurat object to sce object
# sce@colData$orig.ident <- PTtrajectory_A@meta.data$tissue
# sce@colData$orig.ident <- factor(sce@colData$orig.ident, levels = c("Normal","Tumor"))
# # create new meta.data column exp.time and sort levels
# geneFilter <- apply(assays(sce)$counts,1,function(x){
#   sum(x >= 3) >= 10
# })
# sce <- sce[geneFilter, ]
# 
# 
# #  ==== 2.2 normalization ====
# FQnorm <- function(counts){
#   rk <- apply(counts,2,rank,ties.method='min')
#   counts.sort <- apply(counts,2,sort)
#   refdist <- apply(counts.sort,1,median)
#   norm <- apply(rk,2,function(r){ refdist[r] })
#   rownames(norm) <- rownames(counts)
#   return(norm)
# }
# assays(sce)$norm <- FQnorm(assays(sce)$counts)
# # ==== 2.3 dimensionality reduction ====
# # 2.3a) OPTION 1 (PCA)
# pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
# rd1 <- pca$x[,1:2]
# save(rd1, file = 'rd1')
# 
# # 2.3b) OPTION 2 (diffusion maps with destiny)
# library(destiny, quietly = TRUE)
# dm <- DiffusionMap(t(log1p(assays(sce)$norm)))
# rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
# save(rd2, file = 'rd2')
# 
# # add both dimensionality reductions to the sce object
# reducedDims(sce) <- scepleList(PCA = rd1, DiffMap = rd2)
# # ==== 2.4 Clustering Cells ====
# # 2.4a) Add prior clustering information from Seurat to the sce object
# sce@colData$Seurat_clusters <- as.character(Idents(PTtrajectory_A))
# colData(sce)$Seurat_CLUSTERS <- Idents(PTtrajectory_A)
# cl0 <- Idents(PTtrajectory_A)
# save(cl0, file = 'cl0')
# 
# # 2.4b) OPTION 1 Gaussian mixture modeling
# library(mclust, quietly = TRUE)
# cl1 <- Mclust(rd2)$classification #found unique(cl1) = 5 clusters
# save(cl1, file = 'cl1')
# colData(sce)$GMM <- cl1
# 
# # 2.4c) OPTION 2 k means (CAVE: need to arbitrarily choose n of centers!!!)
# cl2 <- kmeans(rd2, centers = 6)$cluster
# save(cl2, file = 'cl2')
# 
# #  ================ 3 USING SLINGSHOT ================
# # ==== 3.1 RUN 1 ====
# # Run Slingshot using the wrapper function "slingshot" (combining "getLineages" and "getCurves")
# # with dimensionality reduction 'DiffMap' and cluster labels identified by 'GMM'.
# sce <- slingshot(sce, reducedDim = 'DiffMap', clusterLabels = 'GMM')
# summary(sce$slingPseudotime_1)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
# 
# pdf('slingshot_DiffMap_GMM.pdf')
# plot(reducedDims(sce)$DiffMap, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
# lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
# legend("right", legend=unique(sce$GMM), 
#        title="GMM Cluster", 
#        col = unique(brewer.pal(9,'Set1')[sce$GMM]), 
#        pch = 20, cex = 1, y.inters=0.7, x.inters=0.5)
# dev.off()
DimPlot(Endo,cols = fibcol)
sce <- as.SingleCellExperiment(Endo)
# Assuming that 'umap' is the key for UMAP in the Seurat object
umap <- Embeddings(Endo, reduction = "umap")
colnames(umap) <- c("UMAP1", "UMAP2")
reducedDim(sce, "UMAP") <- umap
# Assuming 'cell_type_fine' is a column in the metadata of the Seurat object
sce$cell_type_fine <- Endo@meta.data$cell_type_fine

# sim <- slingshot(sce, clusterLabels = 'cell_type_fine', reducedDim = 'UMAP') 
# pdf("trace.pdf",height = 10, width = 10)
# col = brewer.pal(8,'Set1')
# names(col)  = unique(sim$cell_type_fine)
# plot(reducedDims(sim)$UMAP, col = col[sim$cell_type_fine], pch=16, asp = 1) 
# lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
# dev.off()

sim <- slingshot(sce, clusterLabels = 'cell_type_fine', reducedDim = 'UMAP',start.clus = "Fibro-PI16",extend = 'n', stretch = 0)
pdf("trace.pdf",height = 10, width = 10)
col = brewer.pal(8,'Set1')
names(col)  = unique(sim$cell_type_fine)
plot(reducedDims(sim)$UMAP, col = col[sim$cell_type_fine], pch=16, asp = 1) 
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
dev.off()
return(sce)
[1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"



#  ================ 5 IDENTIFY TEMPORALLY EXPRESSED GENES ================
require(gam)
t <- sim$slingPseudotime_1
#  ==== 2.2 normalization ====
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sim)$counts)

# for time, only look at the 100 most variable genes
Y <- log1p(assays(sce)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ to(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#We can then pick out the top genes based on p-value and visualize their expression over developmental time with a heatmap.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sce)$norm[topgenes, order(t, na.last = NA)]
heatclus <- sce$GMM[order(t, na.last = NA)]

pdf('DEG_pseudotime_heatmap.pdf', width=15, height=15)
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus],
        margins = c(5, 5))
dev.off()
gs <- rev(strsplit("ADM, ANGPT2, ANKRD37, APLN, C1QTNF6, CD93, CLDN5, COL4A1, COL4A2, COTL1, CXCR4, DLL4, EDNRB, ESM1, FSCN1, GPIHBP1, HSPG2, IGFBP3, INHBB, ITGA5, JUP, KCNE3, KCNJ8, KDR, LAMA4, LAMB1, LAMC1, LXN, MARCKS, MARCKSL1, MCAM, MEST, MYH9, MYO1B, N4BP3, NID2, NOTCH4, PDGFB, PGF, PLOD1, PLXND1, PMEPA1, PTN, RAMP3, RBP1, RGCC, RHOC, SMAD1, SOX17, SOX4, SPARC, TCF4, UNC5B, VIM",", ")[[1]])
names(gs) <- "tip"
sigs=list("ADM", "ANGPT2", "ANKRD37", "APLN", "C1QTNF6", "CD93", "CLDN5", "COL4A1", "COL4A2", "COTL1", "CXCR4", "DLL4", "EDNRB", "ESM1", "FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", "INHBB", "ITGA5", "JUP", "KCNE3", "KCNJ8", "KDR", "LAMA4", "LAMB1", "LAMC1", "LXN", "MARCKS", "MARCKSL1", "MCAM", "MEST", "MYH9", "MYO1B", "N4BP3", "NID2", "NOTCH4", "PDGFB", "PGF", "PLOD1", "PLXND1", "PMEPA1", "PTN", "RAMP3", "RBP1", "RGCC", "RHOC", "SMAD1", "SOX17", "SOX4", "SPARC", "TCF4", "UNC5B", "VIM")
seu <- AddModuleScore(object = Endo, features = sigs,
                      name = "Exhaustion_Terminal_differentiation", assay = 'RNA', search = F)
FeaturePlot(seu,features = "Exhaustion_Terminal_differentiation1",
            reduction = "umap",pt.size = 1,  cols = c("grey", "red"),combine = FALSE)
pdata <- as.data.frame(seu@meta.data)
C1.meta <- pdata[pdata$Exhaustion_Terminal_differentiation1 == "C1",]
C1.meta <- pdata[pdata$defined == "C1",]
pdata <- as.data.frame(seu@meta.data)
C1.meta <- pdata[pdata$cell_type_fine == "Endo-RGCC",]
C2.meta <- pdata[pdata$cell_type_fine == "Endo-ACKR1",]
C3.meta <- pdata[pdata$cellstate == "CD8-CXCL13",]
C2.meta <- pdata[pdata$cellstate == "CD8-CCL5",]
wilcox.test(C1.meta$Exhaustion_Terminal_differentiation1, C2.meta$Exhaustion_Terminal_differentiation1)
C1.meta
C4.meta <- pdata[pdata$cellstate == "CD8-gdT-TRDC",]
temp_vlnplot <- VlnPlot(object = seu,
                        features = "Exhaustion_Terminal_differentiation1",
                        pt.size = 0,
                        group.by = "cell_type_fine",
                        cols = c("#F7C9CF","#ADD3F4"),
                        ncol = 1)
print(temp_vlnplot)

sigs=list("ACKR1", "AQP1", "C1QTNF9", "CD36", "CSRP2", "EHD4", "FBLN5", "HSPB1", "IL6ST", "JAM2", "LGALS3", "LRG1", "MEOX2", "PLSCR2", "SDPR", "SELP", "SPINT2", "TGFBI", "TGM2", "TMEM176A", "TMEM176B", "TSPAN7", "VWF")
seu <- AddModuleScore(object = Endo, features = sigs,
                      name = "Exhaustion_Terminal_differentiation", assay = 'RNA', search = F)
FeaturePlot(seu,features = "Exhaustion_Terminal_differentiation1",
            reduction = "umap",pt.size = 1,  cols = c("grey", "red"),combine = FALSE)
pdata <- as.data.frame(seu@meta.data)
C1.meta <- pdata[pdata$Exhaustion_Terminal_differentiation1 == "C1",]
C1.meta <- pdata[pdata$defined == "C1",]
pdata <- as.data.frame(seu@meta.data)
C1.meta <- pdata[pdata$cell_type_fine == "Endo-RGCC",]
C2.meta <- pdata[pdata$cell_type_fine == "Endo-ACKR1",]
C3.meta <- pdata[pdata$cellstate == "CD8-CXCL13",]
C2.meta <- pdata[pdata$cellstate == "CD8-CCL5",]
wilcox.test(C1.meta$Exhaustion_Terminal_differentiation1, C2.meta$Exhaustion_Terminal_differentiation1)
C1.meta
C4.meta <- pdata[pdata$cellstate == "CD8-gdT-TRDC",]
temp_vlnplot <- VlnPlot(object = seu,
                        features = "Exhaustion_Terminal_differentiation1",
                        pt.size = 0,
                        group.by = "cell_type_fine",
                        cols = c("#F7C9CF","#ADD3F4") ,
                        ncol = 1)
print(temp_vlnplot)


seurat_10X_integrated <- Endo
temp_exprMatrix <- GetAssayData(object = seurat_10X_integrated,
                                assay = "RNA",
                                slot = "data")
library(Matrix)
library(AUCell)
temp_exprMatrix <-
  Matrix(temp_exprMatrix,
         sparse = T)

# AUCELL
dim(temp_exprMatrix)

print("Buildrankings")
temp_cells_rankings <-
  AUCell_buildRankings(temp_exprMatrix,
                       nCores = 1,
                       plotStats = F)
names(sigs) <- "tip"
temp_cells_AUC <-
  AUCell_calcAUC(geneSets = sigs,
                 rankings = temp_cells_rankings,
                 aucMaxRank = ceiling(0.05 * nrow(temp_cells_rankings)),
                 nCores = 1,
                 verbose = T)
#加载R包library(ggpubr)library(ggimage)library(ggplot2)library(Seurat)
my_comparisons <- list(c("Endo-ACKR1", "Endo-RGCC"))
VlnPlot(Endo, features = "DLL4", group.by = "cell_type_fine")&  theme_bw()&  theme(axis.title.x = element_blank(),        axis.text.x = element_text(color = 'black',face = "bold", size = 12),        axis.text.y = element_text(color = 'black', face = "bold"),        axis.title.y = element_text(color = 'black', face = "bold", size = 15),        panel.grid.major = element_blank(),        panel.grid.minor = element_blank(),        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),        panel.spacing = unit(0.12, "cm"),        plot.title = element_text(hjust = 0.5, face = "bold.italic"),        legend.position = 'none')&  stat_compare_means(method="t.test",hide.ns = F,                      comparisons = my_comparisons,                     label="p.signif",                     bracket.size=0.8,                     tip.length=0,                     size=6)&  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))&   scale_fill_manual(values =  c("#F7C9CF","#ADD3F4"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
A <- singlecell_gene_test(Endo, 
                          genes.use = c('DLL4'),
                          group.by = 'cell_type_fine', 
                          comp = c("Endo-ACKR1", "Endo-RGCC"))
A <- singlecell_gene_test(Endo, 
                          genes.use = c('KDR'),
                          group.by = 'cell_type_fine', 
                          comp = c("Endo-ACKR1", "Endo-RGCC"))

anno_sig <- A$sig
singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = wilcox.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}

plots_violins <- VlnPlot(Endo, 
                         cols = c("#F7C9CF","#ADD3F4"),
                         pt.size = 0,
                         group.by = "cell_type_fine",
                         features = c('KDR'), 
                         ncol = 1, 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Female","Male"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)
anno_sig





VlnPlot(Endo, features = "ESM1", group.by = "cell_type_fine")&  theme_bw()&  theme(axis.title.x = element_blank(),        axis.text.x = element_text(color = 'black',face = "bold", size = 12),        axis.text.y = element_text(color = 'black', face = "bold"),        axis.title.y = element_text(color = 'black', face = "bold", size = 15),        panel.grid.major = element_blank(),        panel.grid.minor = element_blank(),        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),        panel.spacing = unit(0.12, "cm"),        plot.title = element_text(hjust = 0.5, face = "bold.italic"),        legend.position = 'none')&  stat_compare_means(method="t.test",hide.ns = F,                      comparisons = my_comparisons,                     label="p.signif",                     bracket.size=0.8,                     tip.length=0,                     size=6)&  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))&   scale_fill_manual(values = c("#FF5744","#208A42"))



Peri <- subset(int.Mye.RENAMEDIDENT.filt,cell_type_fine == "cDC3-LAMP3")



int.Mye.RENAMEDIDENT.filt$xxx <- Idents(int.Mye.RENAMEDIDENT.filt)
table(int.Mye.RENAMEDIDENT.filt$xxx)
excluded_clusters <- c("T-B", "Fibroblast","Epithelia")
'%notin%' <- Negate('%in%')
int.Mye.RENAMEDIDENT.filt2 <- subset(int.Mye.RENAMEDIDENT.filt, xxx %notin% "excluded_clusters")

