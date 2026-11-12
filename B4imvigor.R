# ============ 加载包 ============
library(IMvigor210CoreBiologies)
library(IOBR)
library(tidyverse)
library(readr)
library(dplyr)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(tibble)
library(tidyr)

# ============ 读取数据 ============
data(cds)
head(counts(cds))
head(fData(cds))
head(pData(cds))

# 提取表达矩阵、注释、临床信息
expMatrix <- counts(cds)
fdata     <- fData(cds)
pdata     <- pData(cds)

# ============ 准备基因注释（entrez_id, length, symbol）============
eff_length2 <- fData(cds)[, c("entrez_id", "length", "symbol")]
rownames(eff_length2) <- eff_length2$entrez_id

# 表达矩阵与注释对齐
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2), ]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

# ============ 计算 TPM ============
x <- expMatrix / eff_length2$length
eset <- t(t(x) / colSums(x)) * 1e6
summary(duplicated(rownames(eset)))

# ============ 基因注释（按symbol去重取均值）============
eset <- IOBR::anno_eset(
  eset       = eset,
  annotation = eff_length2,
  symbol     = "symbol",
  probe      = "entrez_id",
  method     = "mean"
)

# log2 转换
if (max(eset) > 100) eset <- log2(eset + 1)

# ============ 临床信息整理 ============
colnames(pdata) <- gsub(colnames(pdata), pattern = " ", replacement = "_")

pdata <- rownames_to_column(
  pdata[, c("binaryResponse",
            "FMOne_mutation_burden_per_MB",
            "Neoantigen_burden_per_MB",
            "censOS", "os")],
  var = "ID"
)

colnames(pdata) <- c("ID", "BOR_binary", "TMB", "TNB", "status", "time")

# 去除无疗效评价的样本，并转换为 R/NR
pdata <- pdata[!is.na(pdata$BOR_binary), ]
pdata$BOR_binary <- ifelse(pdata$BOR_binary == "CR/PR", "R", "NR")



# CD8-CXCL13 signature
markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")

CD8.CXCL13 <- markers_fine_celltype %>%
  filter(cluster == "CD8-CXCL13" & avg_log2FC >= 2 & p_val_adj < 0.05) %>%
  dplyr::select(gene)

cat("CD8-CXCL13 genes:", nrow(CD8.CXCL13), "\n")

# cDC3-LAMP3 signature
YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")

cdc.lamp3 <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 2 & p_val_adj < 0.05) %>%
  dplyr::select(gene)

cat("cDC3-LAMP3 genes:", nrow(cdc.lamp3), "\n")

# Create gene set list
geneset_list <- list(
  CD8_CXCL13 = CD8.CXCL13$gene,
  cDC_LAMP3 = cdc.lamp3$gene
)


# ============ 2. 封装函数：GSVA + survival +barplot分析 ============

gsvapar <- gsvaParam(eset, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)
library(data.table)

pdat_all <- gsva_es %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%      # 行名变成 ID 列，与 pdata 对应
  left_join(pdata, by = "ID")
pdat_all <- pdat_all %>% filter(!is.na(BOR_binary))

res.cut <- surv_cutpoint(pdat_all, time = "OS.time", 
                         event = "OS_STATUS", 
                         variables = "CD8_CXCL13")
res.cat <- surv_categorize(res.cut)
pdat_all$CD8CXCL13_group<- res.cat[,"CD8_CXCL13"] 
pdat_all$CD8CXCL13_group = factor(pdat_all$CD8CXCL13_group, levels = c("low","high"))

fit = survfit(Surv(OS.time,OS_STATUS) ~ CD8CXCL13_group, data=pdat_all)
p = ggsurvplot(fit, pdat_all, size=0.3, vlegend.labs=unique(pdat_all$CD8CXCL13_group),
               surv.median.line="none", pval=T, conf.int=F,
               palette=c("grey","#00726A"))

### (1) 各组人数、中位生存时间、log-rank P —— surv_summary 一行搞定 ###
surv_summary(fit, data = pdat_all)
# 输出表格中包含：n（各组起始人数）、median（中位生存时间）等

### (2) HR 和 95%CI —— Cox回归 ###
cox_res <- coxph(Surv(OS.time, OS_STATUS) ~ CD8CXCL13_group, data = pdat_all)
summary(cox_res)
# 其中 coef = log(HR)，HR = exp(coef)
# high vs low 的 HR 以 low 组为参照（因为factor levels设了low在前）

# 只提取干净的HR表格
hr_table <- summary(cox_res)$conf.int   # HR、95%CI上下限
hr_table

### (3) log-rank P 值（显式计算，可写入论文）###
sd <- survdiff(Surv(OS.time, OS_STATUS) ~ CD8CXCL13_group, data = pdat_all)
p_value <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_value


# 组装标注文本
hr <- round(hr_table[1, 1], 2)
ci_low <- round(hr_table[1, 3], 2)
ci_high <- round(hr_table[1, 4], 2)
n_low <- table(pdat_all$CD8CXCL13_group)["low"]
n_high <- table(pdat_all$CD8CXCL13_group)["high"]

label_text <- paste0("HR = ", hr, " (", ci_low, "-", ci_high, ")\n",
                     "n(low) = ", n_low, ", n(high) = ", n_high, "\n",
                     "Log-rank P = ", signif(p_value, 3))
label_text
p$plot + 
  annotate("text", x = max(pdat_all$OS.time) * 0.6, 
           y = 0.9, label = label_text, hjust = 0, size = 3.5)

tiff("~/202311_tnbc_caf_05_20240521/GDC.download/KM_IMvigor.tiff",
     width = 120, height = 110, units = "mm", res = 600,   # 单栏或1.5栏尺寸
     compression = "lzw")
print(p)
dev.off()



pre_long <- pdat_all %>%
    select(ID, BOR_binary, CD8_CXCL13, cDC_LAMP3) %>%
    pivot_longer(
        cols = c(CD8_CXCL13, cDC_LAMP3),
        names_to = "signature",
        values_to = "score"
    ) %>%
    mutate(
        signature = factor(signature, 
                           levels = c("CD8_CXCL13", "cDC_LAMP3"),
                           labels = c("CD8-CXCL13", "cDC-LAMP3"))
    )
# Combined plot with facets
p_combined <- ggplot(pre_long, aes(x = BOR_binary, y = score, color = BOR_binary)) +
  geom_boxplot(fill = NA,outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 3.5
  ) +
  facet_wrap(~signature, scales = "free_y") +
  labs(
    title = "IMvigor",
    x = "Response (0=NR, 1=R)",
    y = "GSVA Enrichment Score"
  ) +
  scale_color_manual(values = c("N" = "#737d44", "R" = "#8e2c56")) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  )

print(p_combined)





####survival below###
# Get phenotype data
selected_phenos <- `Melanoma-PRJEB23709.Response (1)`

# Create pdat with all samples
pdat_all <- gsva_es %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample_id') %>%
  left_join(selected_phenos)