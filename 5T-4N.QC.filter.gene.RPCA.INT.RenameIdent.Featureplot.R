gc()
setwd("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812")
CID4495 <- readRDS("~/202311_tnbc_caf_05_20240521/CID4495.RDS")
CID4515 <- readRDS("~/202311_tnbc_caf_05_20240521/CID4515.RDS")
CID44041 <- readRDS("~/202311_tnbc_caf_05_20240521/CID44041.RDS")
CID44991 <- readRDS("~/202311_tnbc_caf_05_20240521/CID44991.RDS")
CID44971 <- readRDS("~/202311_tnbc_caf_05_20240521/CID44971.RDS")
load("~/202311_tnbc_caf_05_20240521/12normal.Object.RData")


#########Qc##########
#########Qc-T##########
samples = read.table("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/5T.txt",header = TRUE)#, stringsAsFactors = FALSE,check.names = FALSE)
minFeature <- 200
maxFeature <- 7500
minCount <- 400
maxCount <- 40000
maxMT <- 10
maxRP <- 50

for (sample_s in samples$sample) {
  print(sample_s)
  sample.tmp.seurat <- get(sample_s)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
                                nCount_RNA < maxCount & percent.mito < maxMT & percent.rp < maxRP & predicted_doublets == F)
  assign(sample_s,sample.tmp.seurat)}

for (sample_s in samples$sample) {
  print(sample_s)
  sample.tmp.seurat <- get(sample_s)
  avg.gene <- median(sample.tmp.seurat$nFeature_RNA)
  print(avg.gene)
}


#########Qc-N##########
samples = read.table("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/4N.txt",header = TRUE)#, stringsAsFactors = FALSE,check.names = FALSE)
minFeature <- 200
maxFeature <- 7500
minCount <- 400
maxCount <- 40000
maxMT <- 10
maxRP <- 50

for (sample_s in samples$sample) {
  print(sample_s)
  sample.tmp.seurat <- get(sample_s)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
                                nCount_RNA < maxCount & percent.mito < maxMT & percent.rp < maxRP)
  assign(sample_s,sample.tmp.seurat)}

for (sample_s in samples$sample) {
  print(sample_s)
  sample.tmp.seurat <- get(sample_s)
  avg.gene <- median(sample.tmp.seurat$nFeature_RNA)
  print(avg.gene)
}


###dim info#########
samples = read.table("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/5T-4N.txt",header = TRUE)#, stringsAsFactors = FALSE,check.names = FALSE)
nCoV.list <- list()
# genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT",
#                   rownames(seurat_object),
#                   value=TRUE, invert=TRUE) #get list of non-ribosomal genes 


for (sample_s in samples$sample) {
  print(sample_s)
  seu <- get(sample_s)
  mito.genes <- grep(pattern = "^MT-", x = rownames(seu), value = TRUE)
  rbl.genes <- grep(pattern = "^RB-", x = rownames(seu), value = TRUE)
  rsl.genes <- grep(pattern = "^RS-", x = rownames(seu), value = TRUE)
  rpl.genes <- grep(pattern = "^RPL-", x = rownames(seu), value = TRUE)
  rbl.genes <- grep(pattern = "^RBL-", x = rownames(seu), value = TRUE)
  rps.genes <- grep(pattern = "^RPS-", x = rownames(seu), value = TRUE)
  rbs.genes <- grep(pattern = "^RBS-", x = rownames(seu), value = TRUE)
  rbl1.genes <- grep(pattern = "^RB", x = rownames(seu), value = TRUE)
  rsl1.genes <- grep(pattern = "^RS", x = rownames(seu), value = TRUE)
  rpl1.genes <- grep(pattern = "^RPL", x = rownames(seu), value = TRUE)
  rbl1.genes <- grep(pattern = "^RBL", x = rownames(seu), value = TRUE)
  rps1.genes <- grep(pattern = "^RPS", x = rownames(seu), value = TRUE)
  rbs1.genes <- grep(pattern = "^RBS", x = rownames(seu), value = TRUE)
  
  counts <- GetAssayData(seu, assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
  assign(sample_s, subset(seu, features = rownames(counts)))
}

for (sample_s in samples$sample) {
  print(sample_s)
  sample.tmp.seurat <- get(sample_s)
  nCoV.list[sample_s] <- sample.tmp.seurat}

result <- do.call(rbind, lapply(nCoV.list, dim))
result
result <- as.data.frame(result)
result
write.table(result,"dim.rbind.5T.4N.xlsx",sep = "\t",quote = F)
#####
setwd("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812")
saveRDS(nCoV.list,"nCoV.list.rds")

s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes

nCoV.list <- SplitObject(seu,split.by="TIG")
nCoV.list <- lapply(X = nCoV.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = nCoV.list,nfeatures = 1000) #Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance

features <- setdiff(features,c(s.genes,g2m.genes))  ##  Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance
# install.packages("future_eapply")  future_lapply
nCoV.list <- lapply(X = nCoV.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = nCoV.list, anchor.features = features, reduction = "rpca")
nCoV.integrated <- IntegrateData(anchorset = immune.anchors)

saveRDS(nCoV.integrated,"nCoV.integrated.rds")

seu <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/seu.rds")
seu[["tissue"]] <- ifelse(seu$TIG %in% c("normal53","normal57","normal62","normal70"),"Normal","Tumor")

###first generate data and scale data in RNA assay
DefaultAssay(seu) <- "RNA"
seu[['percent.mito']] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(seu) <- "integrated"
seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
seu <- RunPCA(seu, verbose = FALSE,npcs = 100)
seu <- ProjectDim(object = seu)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = seu,ndims = 100)
dev.off()
seu = FindNeighbors(seu, assay =  "integrated", dims = 1:20)
seu = FindClusters(object = seu, resolution = c(1, 0.8, 0.6,0.5,0.4,0.3,0.2, 0.1, 0.75, 0.05, 0.023), verbose = TRUE,algorithm = 1, method = "matrix")
seu = RunUMAP(seu, spread = 1.85, dims = 1:15, n.neighbors = 30, min.dist = 0.1)
DimPlot(object = seu, reduction = 'umap',label = FALSE,cols=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#F69B9A",
                                                             "#D65C06","#CC78AA"))
DimPlot(object = data.seurat.filtered, reduction = 'umap',label = TRUE,label.size = 8,repel = T,cols=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#F69B9A",
                                                                                                       "#D65C06","#CC78AA"))
DimPlot(object = data.seurat.filtered, reduction = 'umap',label = TRUE,label.size = 8,repel = T,split.by = "tissue",cols=c("#D9585A","#FEBE80","#B59830","#F38D28","#F69B9A","#4F7CA6","#9FCBE6","#8DCA7D","#F69B9A",
                                                              "#D65C06","#CC78AA"))
data.seurat.filtered <- RenameIdents(seu.with.cluster,
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
markers <- FindAllMarkers(data.seurat.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
write.table(markers,file="markers",sep="\t",quote=F)
library(xlsx,lib.loc = lib_loc)
write.xlsx(markers,file ="sample.marker.xlsx",
           row.names = T,append=T)


















