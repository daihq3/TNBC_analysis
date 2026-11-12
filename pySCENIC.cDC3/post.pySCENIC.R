reg <- read_csv("Seurat_csv/reg.csv")

colnames(reg) <- reg[1,]

reg <- reg[-1,]
reg <- reg[-1,]
colnames(reg)[1] <- "TF"
colnames(reg)[2] <- "MotifID"
# 筛选出 TargetGenes 中包含 "IL15" 的行
# fixed = TRUE 表示匹配固定字符串，速度更快且不涉及正则特殊字符
result <- reg[grepl("IL15", reg$TargetGenes, fixed = TRUE), ]

# 查看结果
print(head(result))

# 如果想查看匹配到了多少行
print(nrow(result))

#97
unique(result$TF)
 [1] "IRF7"   "STAT1"  "ETV7"   "FOXP1"  "MAX"    "ELF1"   "RELB"   "ZEB1"   "ATF4"   "KLF3"   "KLF6"   "KMT2B"  "ZBTB7A" "ZXDC"  
YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")
cdc.lamp3 <- YS.combined_markers %>%
filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05)
result <- reg[grepl("IL15", reg$TargetGenes, fixed = TRUE), ]
write.csv(YS.combined_markers,'dedf.csv')

auc.mtx <- read_csv("Seurat_csv/SCENIC.csv")
# This removes the column AND sets it as rownames in one atomic step
auc.mtx <- column_to_rownames(auc.mtx, var = "Cell")
meta_data <- read_csv("Seurat_csv/meta.data.csv")
regulonActivity_byCellType <- sapply(
    split(meta_data$...1, meta_data$sub_type), # 1. 按 sub_type 分组细胞
    function(cells) {                               # 2. 定义函数：输入是 cells (细胞ID列表)
        # 提取这些细胞在目标 Regulon 中的 AUC 值，并计算平均值
        colMeans(auc.mtx[cells, ], na.rm = TRUE)
    }
)
write.csv(regulonActivity_byCellType,'regulonActivity_byCellType.csv')

