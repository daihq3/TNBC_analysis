library(dplyr)
library(Seurat)

seu <- readRDS("~/202241028/seu.RDS")

# Subset to epithelial cells
# unique(seu@active.ident)
# seu$active.ident <- seu@active.ident
seu$active.ident <- seu$integrated_snn_res.1
DimPlot(seu)
epi <- subset(seu, active.ident %in% c("17","10","5","8","11","13"))
epi <- subset(seu, active.ident %in% c(17,10,5,8,11,13))
FeaturePlot(epi,features = c("EPCAM", "ALDH1A3",  "KRT14","ITGA6","MKI67"))
# Cell cycle assignment
DefaultAssay(epi) <- 'RNA'
epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(epi) <- 'integrated'
epi$cell_cycle <- ifelse(epi$G2M.Score > 0.1 | epi$S.Score > 0.15, 'cycling', 'non-cycling')
DimPlot(epi,group.by = "cell_cycle")

# Fine cell type assignment after manual annotation based on DGE
epi$cell_type_fine <- ifelse(epi$integrated_snn_res.1 %in% c(17,10,5,8,11), 'Luminal','NA') 
epi$cell_type_fine <- ifelse(epi$integrated_snn_res.1 %in% c(13), 'Basal-myoepithelial', epi$cell_type_fine)
epi$barcode <- colnames(epi)
saveRDS(epi, 'epi.cell_type_fine.rds')


# Save object and labels
saveRDS(epi, 'epithelial.rds')
celltype <- epi@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'epithelial_celltype.csv')
