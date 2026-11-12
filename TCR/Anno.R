Combined_MFLN9_20 <- readRDS("~/LAMbrechts_anti_PD1_chemo/0913tnbc/0913doublet/has.integrated.havePCA.haveUMApandTSNE.rds")
# DefaultAssay(CombinedMFLN920) <- 'integrated'
# Combined_MFLN9_20 <- FindClusters(CombinedMFLN920,resolution = c(1.6))
# Combined_MFLN9_20 <- RunUMAP(Combined_MFLN9_20, dims = 1:30) cx
# Idents(CombinedMFLN920) <-  "integrated_snn_res.1.5"
# DimPlot(CombinedMFLN920)
# # markers <- FindAllMarkers(CombinedMFLN920, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = 'MAST',
# #                           assay = 'RNA')
# # write.table(markers,"res1.5.36.cluster.txt")

Idents(Combined_MFLN9_20) <-  "integrated_snn_res.1.5"
table(Idents(Combined_MFLN9_20))


Combined_MFLN9_20 <- RenameIdents(Combined_MFLN9_20, 
                                `1`  = "T/NK",
                                `4`  = "T/NK",
                                `9`  = "T/NK",
                                `25` = "T/NK",
                                `8`  = "T/NK",
                                `0`  = "T/NK",
                                `6`  = "T/NK",
                                `10`  = "T/NK",
                                `12` = "T/NK",
                                `7` = "T/NK",
                                `28` = "T/NK",
                                `22` = "T/NK",
                                `35` = "T/NK",
                                `14` = "T/NK",
                                `31` = "Tumor cell",
                                `23` = "Tumor cell",
                                `5`  = "Tumor cell",
                                `3`  = "Tumor cell",
                                `16` = "Tumor cell",
                                `30` = "Tumor cell",
                                `19` = "Tumor cell", 
                                `21` = "CAF",
                                `18` = "CAF",
                                `26` = "Pericyte",
                                `29` = "B/plasma",
                                `34` = "B/plasma",
                                `2`  = "B/plasma",
                                `27` = "B/plasma",
                                `15` = "B/plasma",
                                `32` = "Myeloid",  
                                `33` = "Myeloid",
                                `17` = "Endothelium",
                                `20` = "Myeloid",
                                `24` = "Myeloid",
                                `11` = "Myeloid",
                                `13` = "Myeloid")

CombinedMFLN920 <- Combined_MFLN9_20
CombinedMFLN920$cell_type_main <- Idents(CombinedMFLN920)
new_levels_order <- c("Tumor cell","Endothelium","CAF","Pericyte","Myeloid","T/NK","B/plasma")
CombinedMFLN920$cell_type_main <- factor(CombinedMFLN920$cell_type_main, levels = new_levels_order)
column.colors <- c(
  "T/NK" = "#D9585A",
  "Tumor cell" = "#FEBE80",
  "CAF" = "#B59830",
  "Myeloid" = "#F38D28",
  "B/plasma" = "#F69B9A",
  "Endothelium" = "#4F7CA6",
  "Pericyte" = "#9FCBE6")

DimPlot(CombinedMFLN920,group.by = "cell_type_main",cols=column.colors,label = T,repel = T)+NoLegend()
# FeaturePlot(CombinedMFLN920,features = "CXCR5")

Idents(CombinedMFLN920) <-  "integrated_snn_res.1.5"
CombinedMFLN920 <- RenameIdents(CombinedMFLN920,
                                `1`  = "CD8",
                                `4`  = "CD8",
                                `9`  = "CD8",
                                `25` = "CD8",
                                `8`  = "CD8",
                                `0`  = "CD4",
                                `6`  = "CD4",
                                `10`  = "CD4",
                                `12` = "CD4",
                                `7` = "CD4",
                                `28` = "CD4",
                                `22` = "Cycling T cell",
                                `35` = "Cycling T cell",
                                `14` = "NK",
                                `31` = "Tumor cell",
                                `23` = "Tumor cell",
                                `5`  = "Tumor cell",
                                `3`  = "Tumor cell",
                                `16` = "Tumor cell",
                                `30` = "Tumor cell",
                                `19` = "Tumor cell",
                                `21` = "CAF",
                                `18` = "CAF",
                                `26` = "Pericyte",
                                `29` = "B cell",
                                `34` = "B cell",
                                `2`  = "B cell",
                                `27` = "Plasma cell",
                                `15` = "Plasma cell",
                                `32` = "pDC",
                                `33` = "Mast cell",
                                `17` = "Endothelium",
                                `20` = "TAM",
                                `24` = "mDC",
                                `11` = "TAM",
                                `13` = "Monocyte")



table(Idents(CombinedMFLN920))
length(unique(Idents(CombinedMFLN920)))
#23
# CD8-GZMK     CD8-CXCL13       CD8-IL7R      CD8-ISG15   CD8-gdT-TRDC       CD4-CCR7 
# 3706           2382           1693            634           2016           5046 
# CD4-FOXP3     CD4-CXCL13      CD4-ISG15 Cycling T cell             NK     Tumor cell 
# 3933           3534            511            929           1270           8390 
# iCAF          myCAF       Pericyte         B cell    Plasma cell            pDC 
# 854            923            582           3322           1773            225 
# Mast cell    Endothelium     Macrophage            mDC       Monocyte 
# 128            931           2505            676           1465 
CombinedMFLN920$cell_type_intermediate <- Idents(CombinedMFLN920) 
unique(CombinedMFLN920$cell_type_intermediate)
new_levels_order <- c("Tumor cell","Endothelium","CAF","Pericyte","Mast cell","pDC","mDC","Monocyte","TAM","CD4","CD8","NK","Cycling T cell","B cell","Plasma cell")
CombinedMFLN920$cell_type_intermediate <- factor(CombinedMFLN920$cell_type_intermediate, levels = new_levels_order)
DimPlot(CombinedMFLN920,group.by = "cell_type_intermediate")

Idents(CombinedMFLN920) <-  "integrated_snn_res.1.5"
CombinedMFLN920 <- RenameIdents(CombinedMFLN920,
                                `1`  = "CD8-GZMK",
                                `4`  = "CD8-CXCL13",
                                `9`  = "CD8-IL7R",
                                `25` = "CD8-ISG15",
                                `8`  = "CD8-gdT-TRDC",
                                `0`  = "CD4-CCR7",
                                `6`  = "CD4-FOXP3",
                                `10`  = "CD4-FOXP3",
                                `12` = "CD4-CXCL13",
                                `7` = "CD4-CXCL13",
                                `28` = "CD4-ISG15",
                                `22` = "Cycling T cell",
                                `35` = "Cycling T cell",
                                `14` = "NK",
                                `31` = "Tumor cell",
                                `23` = "Tumor cell",
                                `5`  = "Tumor cell",
                                `3`  = "Tumor cell",
                                `16` = "Tumor cell",
                                `30` = "Tumor cell",
                                `19` = "Tumor cell",
                                `21` = "iCAF",
                                `18` = "myCAF",
                                `26` = "Pericyte",
                                `29` = "B cell",
                                `34` = "B cell",
                                `2`  = "B cell",
                                `27` = "Plasma cell",
                                `15` = "Plasma cell",
                                `32` = "pDC",
                                `33` = "Mast cell",
                                `17` = "Endothelium",
                                `20` = "Macrophage",
                                `24` = "mDC",
                                `11` = "Macrophage",
                                `13` = "Monocyte")

CombinedMFLN920$cell_type_fine <- Idents(CombinedMFLN920) 
unique(CombinedMFLN920$cell_type_fine)
Mye <- subset(CombinedMFLN920, cell_type_fine %in% c("pDC","DC",
                                                     "Mast cell",
                                                     "Macrophage",
                                                     "mDC",
                                                     "Macrophage",
                                                     "Monocyte"))

# new_levels_order <- c("Tumor cell","Endothelium","CAF","Pericyte","Mast cell","pDC","mDC","Monocyte","TAM","CD4","CD8","NK","Cycling T cell","B cell","Plasma cell")
# CombinedMFLN920$cell_type_intermediate <- factor(CombinedMFLN920$cell_type_intermediate, levels = new_levels_order)
# DimPlot(CombinedMFLN920,group.by = "cell_type_intermediate")
# saveRDS(seu,"~/202241028/Lambrecht.TNBC/Combined_TNBC_Lambrecht.rescale.rds")
`CombinedMFLN920` <- readRDS("/home/data/t180330/202241028/Lambrecht.TNBC.02/Mye/Combined_TNBC_Lambrecht.20241109.rds")
CombinedMFLN920 <- Combined_TNBC_Lambrecht
Idents(CombinedMFLN920) <- "cell_type_main"
levels(
  Idents(CombinedMFLN920) 
)
markers <- FindAllMarkers(CombinedMFLN920, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = 'MAST',
                          assay = 'RNA')

markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) -> top10
write.table(top10,"top30_markers_cell_type_main_Lambrecht.rescale.tsv",sep = '\t')
write.table(markers,"markers_cell_type_main_Lambrecht.rescale.tsv",sep = '\t')
Idents(CombinedMFLN920) <- "cell_type_intermediate"
levels(
  Idents(CombinedMFLN920) 
)
markers <- FindAllMarkers(CombinedMFLN920, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = 'MAST',
                          assay = 'RNA')
markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) -> top10
write.table(top10,"top30_markers_cell_type_intermediate_Lambrecht.rescale.tsv",sep = '\t')
write.table(markers,"markers_cell_type_intermediate_Lambrecht.rescale.tsv",sep = '\t')


Idents(CombinedMFLN920) <- "cell_type_fine"
levels(
  Idents(CombinedMFLN920) 
)
options(future.globals.maxSize = Inf)
DefaultAssay(CombinedMFLN920) <- "RNA"
markers <- FindAllMarkers(CombinedMFLN920, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0,test.use = 'MAST',
                          assay = 'RNA')
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
write.table(top10,"top10_markers_cell_type_fine_Lambrecht.rescale.tsv",sep = '\t')
write.table(markers,"markers_cell_type_fine_Lambrecht.rescale.tsv",sep = '\t')
cell_colors <- cluter_mye_cols
cluster.averages_Pan02_CD8 <- AverageExpression(CombinedMFLN920, return.seurat = TRUE)
ph <- DoHeatmap(cluster.averages_Pan02_CD8, 
                features = top10$gene, 
                group.colors = column.colors3,
                hjust = F,
                draw.lines = F,
                angle = 45,
                size = 4) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) +
  theme(axis.text.y = element_blank())

output_file <- "~/202241028/Lambrecht.TNBC.02/Mye/03/heatmap.pdf"
# 使用pdf函数打开一个PDF图形设备
pdf(file = output_file, width = 10, height = 10)
# 显示热图
plot(ph)
# 关闭图形设备
dev.off()
nCoV_markers_rev <- c(
"CLEC9A", "FLT3", "IDO1",
"CD1C", "FCER1A", "HLA-DQA1",
"LAMP3", "CCR7", "BIRC3",
"LILRA4", "JCHAIN", "GZMB",
"FCN1", "S100A8", "S100A9",
"SPP1", "CTSL", "MARCO",
"STAB1", "FOLR2", "SELENOP",
"IL1B", "NLRP3", "CLEC5A",
"C1QC", "C1QA", "APOE",
"CXCL10", "CXCL11", "ISG15",
"TPSAB1", "KIT", "CPA3",
"MKI67", "CDC20", "TOP2A")
Idents(harmonized_seurat2) = harmonized_seurat2$cell_type_fine
library(dplyr)
library(forcats)
major_celltype_order <- rev(levels(harmonized_seurat2$cell_type_fine))
harmonized_seurat2@meta.data$cell_type_fine <- fct_relevel(harmonized_seurat2@meta.data$cell_type_fine, major_celltype_order)
pdf(file="tcr-heatmap-seurat.pdf", width = 9, height = 3.5)
pp = DotPlot(harmonized_seurat2, features = nCoV_markers_rev,cols = c('white','#F8766D'),dot.scale=4.5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 1))
print(pp)
dev.off()

p1 <- pp
exp <- p1$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)
library(ggplot2)


harmonized_seurat2 <- RenameIdents(YS.combined,
                                   "15"="cDC1-CLEC9A",     # DC细胞群
                                   "6"="cDC2-CD1C",        # DC细胞群
                                   "13"="cDC-LAMP3",      # DC细胞群
                                   "8"="pDC-LILRA4",      # DC细胞群
                                   "9"="Mono-FCN1",       # 单核细胞
                                   "11"="Macro-SPP1",      # 巨噬细胞
                                   "3"="Macro-SPP1",      # 巨噬细胞
                                   "5"="Macro-SPP1",      # 巨噬细胞
                                   "4"="Macro-FOLR2",      # 巨噬细胞
                                   "1"="Macro-IL1B",     # 巨噬细胞
                                   "2"="Macro-C1QC",     # 巨噬细胞
                                   "0"="Macro-C1QC",    # 巨噬细胞
                                   "7"="Macro-CXCL10",                                  
                                   "10"="Mast-TPSAB1", 
                                   "14"="Cycling Myeloid", 
                                   "12"="remove")   # 巨噬细胞
# list_genes <- list(
#   "cDC1-CLEC9A" = c("CLEC9A", "FLT3", "IDO1"),
#   "cDC2-CD1C" = c("CD1C", "FCER1A", "HLA-DQB1"), # 修正了等号
#   "cDC-LAMP3" = c("LAMP3", "CCR7", "BIRC3"),
#   "pDC-LILRA4" = c("LILRA4", "JCHAIN", "GZMB"),
#   "Mono-FCN1" = c("FCN1", "S100A8", "S100A9"),
#   "Macro-SPP1" = c("SPP1", "CTSB", "MACRO"),
#   "Macro-FOLR2" = c("LYVE1", "FOLR2", "SELENOP"), # 修正了基因名
#   "Macro-IL1B" = c("IL1B", "NLRP3", "EREG"),
#   "Macro-C1QC" = c("C1QC", "C1QA", "APOE"),
#   "Macro-CXCL10" = c("CXCL10", "CXCL11", "ISG15"),
#   "Mast-TPSAB1" = c("TPSAB1", "KIT", "CPA3"),
#   "Cycling Myeloid" = c("MKI67", "CDC20", "TOP2A") # 修正了基因名
# )


library(RColorBrewer)
# #----------------------------------------------------------------------------------------------------------
# #### Fig.S1A ####
# #----------------------------------------------------------------------------------------------------------
# MHC2 <- c("CD74","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1",
#           "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2",
#           "HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5")
# MHC2 <- nCoV_markers_rev
# MHC2 <- MHC2[MHC2 %in% rownames(harmonized_seurat2)]
# DotPlot(object = harmonized_seurat2, features = MHC2,scale = T,group.by = "cell_type_fine") + 
#   scale_colour_gradientn(colors=brewer.pal(10, "BuPu")) + theme_bw() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_x_discrete(breaks=MHC2,labels=MHC2)#MHC2

pdf(file="tcr-heatmap-seurat.pdf", width = 12, height = 3.5)
pp = DotPlot(harmonized_seurat2, features = nCoV_markers_rev,cols = c('white', "#308a46"),dot.scale=4.5,group.by = "cell_type_fine") + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 1))
print(pp)
dev.off()
IL1B

# aaaaa -------------------------------------------------------------------
# 
# nonono ------------------------------------------------------------------



# #----------------------------------------------------------------------------------------------------------
# #### Fig.10B ####
# #----------------------------------------------------------------------------------------------------------
# load("./Revision/Github/Input_data/MarkerExp.rda")
# load("./Revision/Github/Input_data/cellType_All.rda")
# df <- unique(Bcell@meta.data[,c("orig.ident","Field")]);rownames(df) <- df$orig.ident
# cellType_All$Field <- df[cellType_All$orig.ident,"Field"]
# obj <- CreateSeuratObject(t(MarkerExp),assay = "RNA",meta.data = cellType_All)
# obj@assays$RNA@data <- as.matrix(t(MarkerExp))
# obj <- subset(obj,subset=cellType %in% c("CD4T","CD8T","Fibroblast","Myeloid"))
# obj$cellType <- paste0(obj$Field,"-",obj$cellType)
# p1 <- DotPlot(object = obj, features = colnames(MarkerExp),assay = "RNA",scale = T,group.by = "cellType") + 
#   scale_colour_gradientn(colors=brewer.pal(9, "BuPu"))  + theme_bw() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_x_discrete(breaks=colnames(MarkerExp),labels=colnames(MarkerExp))





# (harmonized_seurat2$cell_type_fine)
# unique(harmonized_seurat2$cell_type_fine)
DimPlot(harmonized_seurat2, group.by = "cell_type_fine", cols = c("#A6C694", "#4DB7E2", "#FEE2D5", "#86B5D4", "#BFE3ED", "#CAE5D7", "#ACA7C2", "#61C1C5", "plum3", "#fdca9f", "#d6c8e2", "#a9dd9b", "cyan3", "#d2b0ed", "#8f657d", "#9b9349", "#E4556E", "#f1989c", "#7f3d7e", "#34bfd6",
                                                                       "#7f7e81",
                                                                       "#f4ac71",
                                                                       "#4b94c1", "brown", "lightblue", "cadetblue", "cyan4", "chartreuse3", "chocolate", "coral1", "darkorange",
                                                                       "cornflowerblue", "darkgoldenrod", "darkolivegreen"), label = TRUE, repel = TRUE) + NoLegend() ->
  plot_intermediate
ggsave("cell_type_fine_dimplot.pdf", plot_intermediate, width = 10, height = 10)

# 第三个DimPlot，保存为cell_type_main_dimplot.pdf

