library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(GSVA)

exprSet <- data.table::fread("~/202311_tnbc_caf_05_20240521/GDC.download/brca_metabric/data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",data.table = F)

# library(genefu)
# data(pam50.robust) #必要加载步骤

exprSet <- exprSet[,-2]

duplicated(exprSet$Hugo_Symbol) %>% table()  

exprSet <- limma::avereps(exprSet[,-1],ID=exprSet$Hugo_Symbol)  
head(exprSet [,1:6])
exprSet <- as.data.frame(exprSet)  

boxplot(exprSet[,1:6],las=2)
expr <- exprSet
cl <- data.table::fread("~/202311_tnbc_caf_05_20240521/GDC.download/brca_metabric/data_clinical_sample.txt",data.table = F,skip = 4)
OS <- data.table::fread("~/202311_tnbc_caf_05_20240521/GDC.download/brca_metabric/data_clinical_patient.txt")
names(cl)
# 假设您的原始数据框是 df
# df <- read.csv("your_file.csv") 或者其他方式创建的

# 定义新的列名向量
new_colnames <- c(
  "PATIENT_ID", "SAMPLE_ID", "CANCER_TYPE",
  "CANCER_TYPE_DETAILED", "ER_STATUS", "HER2_STATUS",
  "GRADE", "ONCOTREE_CODE", "PR_STATUS",
  "SAMPLE_TYPE", "TUMOR_SIZE", "TUMOR_STAGE",
  "TMB_NONSYNONYMOUS"
)

# 直接将新列名赋值给数据框的列名属性
colnames(cl) <- new_colnames

# 检查一下结果
head(cl)
pd <- cl %>% filter(ER_STATUS!=""&  
                      HER2_STATUS!=""&  
                      PR_STATUS!="")  

# install.packages("dplyr") # 如果未安装，请运行此行
library(dplyr)

# 同样，先确保 GRADE 是数值型
pd <- pd %>%
  mutate(GRADE_NUM = as.numeric(as.character(GRADE)))

# 使用 case_when 进行分类
pd <- pd %>%
  mutate(Type = case_when(
    # TNBC: ER、PR、HER2均阴性
    ER_STATUS == "Negative" & PR_STATUS == "Negative" & HER2_STATUS == "Negative" ~ "TNBC",
    
    # HER2型: ER和PR阴性，HER2阳性
    ER_STATUS == "Negative" & PR_STATUS == "Negative" & HER2_STATUS == "Positive" ~ "HER2",
    
    # Luminal B型: ER阳性、PR阴性、HER2阳性
    ER_STATUS == "Positive" & PR_STATUS == "Negative" & HER2_STATUS == "Positive" ~ "Luminal B",
    
    # Luminal A型: ER或PR阳性，HER2阴性
    (ER_STATUS == "Positive" | PR_STATUS == "Positive") & HER2_STATUS == "Negative" ~ "Luminal A",
    
    # 默认值: 以上条件都不满足
    TRUE ~ "Unknown"
  ))
table(pd$Type)
# 提取 TNBC (对应之前的 choose_samples1)
choose_samples1 <- pd %>%
  filter(Type == "TNBC") %>%
  pull(PATIENT_ID)

# 提取 HER2 (对应之前的 choose_samples2)
choose_samples2 <- pd %>%
  filter(Type == "HER2") %>%
  pull(PATIENT_ID)

# 提取 Luminal B (对应之前的 choose_samples3)
choose_samples3 <- pd %>%
  filter(Type == "Luminal A") %>%
  pull(PATIENT_ID)

# 提取 Luminal A (对应之前的 choose_samples4)
choose_samples4 <- pd %>%
  filter(Type == "Luminal B") %>%
  pull(PATIENT_ID)

# 3. 打印数量进行核对
cat("TNBC samples:", length(choose_samples1), "\n")
cat("HER2 samples:", length(choose_samples2), "\n")
cat("Luminal A samples:", length(choose_samples3), "\n")
cat("Luminal B samples:", length(choose_samples4), "\n")



# CD8-CXCL13 signature
markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")

CD8.CXCL13 <- markers_fine_celltype %>%
  filter(cluster == "CD8-CXCL13" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)

cat("CD8-CXCL13 genes:", nrow(CD8.CXCL13), "\n")

# cDC3-LAMP3 signature
YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")

cdc.lamp3 <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)

cat("cDC3-LAMP3 genes:", nrow(cdc.lamp3), "\n")

# Create gene set list
geneset_list <- list(
  CD8_CXCL13 = CD8.CXCL13$gene,
  cDC_LAMP3 = cdc.lamp3$gene
)



expr <- as.matrix(expr)
gsvapar <- gsvaParam(expr, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)

all_samples <- c(choose_samples1, choose_samples2, choose_samples3, choose_samples4)
sample_types <- c(rep("Basal", length(choose_samples1)), 
                  rep("Her2", length(choose_samples2)),
                  rep("LumA", length(choose_samples3)),
                  rep("LumB", length(choose_samples4)))

sample_type_df <- data.frame(case = all_samples, Types = sample_types)
print(sample_type_df)
pdat <- gsva_es %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('case') %>%
  left_join(sample_type_df, by = "case") %>%
  filter(!is.na(Types)) %>%
  mutate(Types = factor(Types, levels = c("Basal", "Her2", "LumA", "LumB")))
bxp <- ggboxplot(pdat, x = "Types", y = "CD8_CXCL13", fill = "Types",
                 palette = c("#009fdf", "#ffa100", "#8e1e98", "#c68ecb"),
                 ylab = "Signature Score")
print(bxp)

kruskal.test(CD8_CXCL13 ~ Types, data = pdat)
Basal_her2 <- subset(liu_clin_sig.df, Types %in% c("Basal", "Her2"))
wilcox.test(score ~ Types, data = Basal_her2)
Basal_her2 <- subset(liu_clin_sig.df, Types %in% c("Basal", "LumA"))
wilcox.test(score ~ Types, data = Basal_her2)
Basal_her2 <- subset(liu_clin_sig.df, Types %in% c("Basal", "LumB"))
wilcox.test(score ~ Types, data = Basal_her2)

pdat2 <- pdat %>%
    left_join(cl, by = c("case" = "PATIENT_ID")) %>%
    mutate(TUMOR_STAGE = as.numeric(TUMOR_STAGE)) %>%
    filter(TUMOR_STAGE %in% c(1, 2, 3, 4))
pdat2 <- pdat2 %>%
  mutate(TUMOR_STAGE = factor(TUMOR_STAGE, levels = c(1, 2, 3, 4)))
stage_colors <- c("1" = "#a6cee3",
                  "2" = "#1f78b4",
                  "3" = "#b2df8a",
                  "4" = "#33a02c")

bxp <- ggplot(pdat2, aes(x = TUMOR_STAGE, y = CD8_CXCL13, fill = TUMOR_STAGE)) +
    geom_jitter(aes(color = TUMOR_STAGE),
                width = 0.15,
                alpha = 0.6,
                size = 1.5) +
    geom_boxplot(alpha = 0.4,
                 outlier.shape = NA,
                 width = 0.7) +
    scale_fill_manual(values = stage_colors) +
    scale_color_manual(values = stage_colors) +
    labs(title = "All subtypes", x = "Stage", y = "FOLR2") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
bxp <- bxp + stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
        c("1", "2"),
        c("1", "3"),
        c("1", "4")
    ),
    label = "p.format",
    hide.ns = FALSE
)

library(ggplot2)
library(ggpubr)

# 亚型因子顺序 + 颜色（与原图一致）
pdat$Types <- factor(pdat$Types,
                     levels = c("Basal", "Her2", "LumA", "LumB"))
cols <- c("Basal" = "#C7102E",     # 深红
          "Her2" = "#E6399B",      # 品红
          "LumA"  = "#000bc4",     # 深蓝
          "LumB"  = "#2e83b6")  # 浅蓝

# 绘图
p <- ggplot(pdat, aes(x = CD8_CXCL13, y = cDC_LAMP3, color = Types)) +
    geom_point(size = 0.4, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = "Human CD8_CXCL13\nGSVA Enrichment Score",
         y = "GSVA Enrichment Score\nHuman cDC_LAMP3") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",          # <-- 去掉图例
          axis.title = element_text(size = 14, face = "bold"),
          plot.tag = element_text(face = "bold", size = 18)) +
    labs(tag = "H")
p

# 手动计算每个亚型的 Pearson R 和 p，用于图例/标注
stats <- do.call(rbind, lapply(split(pdat, pdat$Types), function(d) {
    ct <- cor.test(d$CD8_CXCL13, d$cDC_LAMP3, method = "pearson")
    data.frame(PAM50 = unique(d$Types),
               R = round(ct$estimate, 2),
               p = ct$p.value)
}))
print(stats)
# 若想显示成图中那种 R_x.xx, p<2.2x10^-16 的形式，可以自己格式化：
# sprintf("R=%.2f, p=%.1g", stats$R, stats$p)
tiff("~/202311_tnbc_caf_05_20240521/GDC.download/cor.all.subtypes.tiff",
     width = 100, height = 100, units = "mm", res = 600,   # 单栏或1.5栏尺寸
     compression = "lzw")
print(p)
dev.off()
