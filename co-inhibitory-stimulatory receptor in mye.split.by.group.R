Mye.filt$celltype <- Idents(Mye.filt)

sel.genes <- c('TNFSF4', 'CD80', 'TNFSF15', 'ICOSLG', 'CD40', 'CD86', 'TNFSF18', 'CD274','PDCD1LG2','LGALS9', 'BTNL2', 'HHLA2', 'CD276', 'VTCN1', 'CD200','CD200R1')


# Define custom gene list
# custom_genes <- c('TNFSF4', 'CD80', 'TNFSF15', 'ICOSLG', 'CD40', 'CD86', 'TNFSF18', 'PDCD1LG2', 'LGALS9', 'BTNL2', 'HHLA2', 'CD276', 'VTCN1', 'CD274', 'CD200')

Mye.filt@meta.data$celltype_condition=paste(
  as.character(Mye.filt@meta.data$celltype),
  Mye.filt@meta.data$tissue,
  sep = "_"
)
celltype_order <- c("cDC1-CLEC9A", "cDC2-CD1C", "cDC3-LAMP3", "pDC-LILRA4", 
                    "Mono-FCN1", "Macro-SPP1", "Macro-C1QC", "Macro-FOLR2", "Macro-CXCL10")
sampletype_order <- c("Normal", "Tumor")
# Create the full order list with "Normal" and "Tumor" appended to each cell type
celltype_condition_order <- c(
  sapply(sampletype_order, function(st) paste(st, celltype_order, sep = "_"))
)

# 假设 Mye.filt@meta.data 已经包含了 celltype 和 tissue 列
# 创建 celltype_condition 列，并转换为因子
celltype_condition_order = c(
  "cDC1-CLEC9A_Normal", "cDC1-CLEC9A_Tumor", "cDC2-CD1C_Normal", "cDC2-CD1C_Tumor",
  "cDC3-LAMP3_Normal", "cDC3-LAMP3_Tumor", "pDC-LILRA4_Normal", "pDC-LILRA4_Tumor",
  "Mono-FCN1_Normal", "Mono-FCN1_Tumor", "Macro-SPP1_Normal", "Macro-SPP1_Tumor",
  "Macro-C1QC_Normal", "Macro-C1QC_Tumor", "Macro-FOLR2_Normal", "Macro-FOLR2_Tumor",
  "Macro-CXCL10_Normal", "Macro-CXCL10_Tumor")


Mye.filt@meta.data$celltype_condition=factor(
  Mye.filt@meta.data$celltype_condition,
  levels = celltype_condition_order)




Idents(Mye.filt) <- "celltype_condition"
p <- DotPlot(object = Mye.filt,features = rev(sel.genes), assay = "RNA", scale = T) +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.text.x = element_text(size = 12)) +  coord_flip() +
  scale_x_discrete(breaks = sel.genes, labels = sel.genes)
p 


sel.genes <- c('CXCL16')
p <- VlnPlot(object = Mye.filt,features = rev(sel.genes),cols = c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"))
Idents(Mye.filt) <- "celltype_condition"
p <- DotPlot(object = Mye.filt,features = rev(sel.genes), assay = "RNA", scale = T) +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.text.x = element_text(size = 12)) +  coord_flip() +
  scale_x_discrete(breaks = sel.genes, labels = sel.genes)
p 

sel.genes <- c("CXCR6","CCR2","CCR5","CXCR3","CCR1","CCR8","CCR3","CCR6","CXCR4","XCR1","CXCR5","CCR9","CCR7","CCR4","CXCR2","CXCR1","CCR10","CX3CR1","CXCR7")
p <- DotPlot(object = Tcell.filt,features = rev(sel.genes), assay = "RNA", scale = T) +
  scale_colour_gradientn(colors = brewer.pal(9, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_text(size = 12))  +
  theme(axis.text.x = element_text(size = 12)) +  coord_flip() +
  scale_x_discrete(breaks = sel.genes, labels = sel.genes)
p 
