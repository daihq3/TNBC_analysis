# Melanoma PRJEB23709 GSVA Analysis with CD8-CXCL13 and cDC-LAMP3 signatures

library(readr)
library(dplyr)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(tibble)
library(tidyr)

cat("=== Melanoma GSVA Analysis ===\n")

# ============================================
# STEP 1: Load gene signatures
# ============================================

cat("\n1. Loading gene signatures...\n")

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

# ============================================
# STEP 2: Process expression data
# ============================================

cat("\n2. Processing expression data...\n")

# Load expression data
gexp.df <- `Melanoma-PRJEB23709.Response`

# Set rownames and transpose
rowne <- gexp.df$GENE_SYMBOL
gexp.df <- gexp.df[, -1]  # Remove GENE_SYMBOL column
rownames(gexp.df) <- rowne

gexp.df <- as.matrix(gexp.df)

# Convert to numeric
gexp.clean <- apply(gexp.df, 2, as.numeric)
rownames(gexp.clean) <- rowne

# Transpose for GSVA (genes as rows, samples as columns)

cat("Expression matrix dimensions:", dim(gexp.matrix), "(samples x genes)\n")
cat("Genes in matrix:", length(intersect(rownames(gexp.clean), c(CD8.CXCL13$gene, cdc.lamp3$gene))), "overlap with signatures\n")

# ============================================
# STEP 3: Run GSVA
# ============================================

cat("\n3. Running GSVA...\n")

gsvapar <- gsvaParam(gexp.clean, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)

cat("GSVA complete. Output dimensions:", dim(gsva_es), "\n")

# ============================================
# STEP 4: Merge with phenotype data
# ============================================

cat("\n4. Merging with phenotype data...\n")

selected_phenos <- `Melanoma-PRJEB23709.Response (1)`

pdat <- gsva_es %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample_id') %>%
  left_join(selected_phenos, by = "sample_id")

cat("Total samples:", nrow(pdat), "\n")
cat("Treatment groups:", table(pdat$Treatment), "\n")
cat("Response groups:", table(pdat$response_NR), "\n")

# ============================================
# STEP 5: Filter PRE treatment only
# ============================================

cat("\n5. Filtering PRE treatment...\n")

pre_pdat <- pdat %>% filter(Treatment == "PRE")

cat("PRE samples:", nrow(pre_pdat), "\n")
cat("PRE response distribution:\n")
print(table(pre_pdat$response_NR))

# ============================================
# STEP 6: Create boxplots
# ============================================

cat("\n6. Creating boxplots...\n")

# Function to create boxplot with stats
create_melanoma_plot <- function(data, signature_col, title, ylabel) {
  data <- data %>%
    mutate(response_NR = factor(response_NR, levels = c("N", "R")))
  
  # Wilcoxon test (non-parametric, better for small samples)
  test_result <- wilcox.test(as.formula(paste(signature_col, "~ response_NR")), data = data)
  pval <- signif(test_result$p.value, 4)
  
  # Get sample sizes
  n_table <- data %>%
    group_by(response_NR) %>%
    summarise(n = n(), .groups = "drop")
  
  cat(title, "- Wilcoxon p-value:", pval, "\n")
  cat("  Non-responders (0): n =", n_table$n[n_table$response_NR == "N"], 
      ", mean =", round(mean(data[[signature_col]][data$response_NR == "R"], na.rm = TRUE), 3), "\n")
  cat("  Responders (1): n =", n_table$n[n_table$response_NR == "R"], 
      ", mean =", round(mean(data[[signature_col]][data$response_NR == "R"], na.rm = TRUE), 3), "\n")
  
  # Create plot
  p <- ggplot(data, aes(x = response_NR, y = .data[[signature_col]], color = response_NR)) +
    geom_boxplot(fill = NA, width = 0.6, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.x.npc = "center",
      label.y.npc = "top",
      size = 4,
      vjust = -0.5
    ) +
    labs(
      title = title,
      x = "Response (0=N, 1=R)",
      y = ylabel
    ) +
    scale_color_manual(values = c("N" = "#737d44", "R" = "#8e2c56")) +
    scale_x_discrete(labels = c("N" = paste0("Non-responder\nn=", n_table$n[n_table$response_NR == "N"]),
                                "R" = paste0("Responder\nn=", n_table$n[n_table$response_NR == "R"]))) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )

  return(p)
}

# CD8-CXCL13 plot
p1 <- create_melanoma_plot(
  pre_pdat, 
  "CD8_CXCL13", 
  "CD8-CXCL13 Signature - Melanoma PRE",
  "GSVA Enrichment Score"
)
print(p1)

# cDC-LAMP3 plot
p2 <- create_melanoma_plot(
  pre_pdat,
  "cDC_LAMP3",
  "cDC-LAMP3 Signature - Melanoma PRE",
  "GSVA Enrichment Score"
)
print(p2)

# ============================================
# STEP 7: Combined plot
# ============================================

cat("\n7. Creating combined plot...\n")

# Reshape data for combined plot
pre_long <- pre_pdat %>%
  select(sample_id, response_NR, CD8_CXCL13, cDC_LAMP3) %>%
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
p_combined <- ggplot(pre_long, aes(x = response_NR, y = score, color = response_NR)) +
  geom_boxplot(fill = NA,outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 3.5
  ) +
  facet_wrap(~signature, scales = "free_y") +
  labs(
    title = "Melanoma PRE Treatment - Biomarker Signatures",
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

cat("\n=== Analysis Complete ===\n")


####survival below###
# Get phenotype data

# Filter for PRE-treatment samples
selected_phenos <- `Melanoma-PRJEB23709.Response (1)`

markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")
selected_genes <- markers_fine_celltype %>%
  filter(cluster == "CD8-CXCL13" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)

library(dplyr)
geneset_list <- list(IFNG_signature = selected_genes$gene)

YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")


selected_genes <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)
#120
library(dplyr)
geneset_list=list(IFNG_signature=selected_genes$gene)


# ============ 1. 先分别提取 PRE 和 EDT 样本 ============
pheno_pre  <- selected_phenos %>% filter(Treatment == "PRE")
pheno_post <- selected_phenos %>% filter(Treatment == "EDT")

# 按样本名提取对应的表达矩阵子集（假设 gexp.clean 的列名是样本ID）
gexp_pre  <- gexp.clean[, colnames(gexp.clean) %in% pheno_pre$sample_id]
gexp_post <- gexp.clean[, colnames(gexp.clean) %in% pheno_post$sample_id]

# ============ 2. 封装函数：GSVA + 生存分析 ============
run_survival_analysis <- function(gexp_sub, pheno_sub, title_prefix) {
  
  if (ncol(gexp_sub) < 10) {
    cat("\n===", title_prefix, "===\n")
    cat("Insufficient samples:", ncol(gexp_sub), "\n")
    return(NULL)
  }
  
  # ---- GSVA 打分（在函数内部进行）----
  gsvapar <- gsvaParam(gexp_sub, geneset_list, maxDiff = TRUE)
  gsva_es <- gsva(gsvapar)
  
  pdat <- gsva_es %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('sample_id') %>%
    left_join(pheno_sub, by = "sample_id")
  
  # ---- 生存数据准备 ----
  pdat <- pdat %>%
    mutate(BR2 = ifelse(`vital status` == "Dead", 1, 0))
  pdat$OS.time <- pdat$`overall survival (days)` / 365
  
  # ---- 最佳切点 ----
  res.cut <- surv_cutpoint(pdat, time = "OS.time",
                           event = "BR2",
                           variables = "IFNG_signature")
  res.cat <- surv_categorize(res.cut)
  
  pdat$CDC3_group <- res.cat[, "IFNG_signature"]
  pdat$CDC3_group <- factor(pdat$CDC3_group, levels = c("low", "high"))
  
  # ---- KM 曲线 ----
  fit <- survfit(Surv(OS.time, BR2) ~ CDC3_group, data = pdat)
  p <- ggsurvplot(fit, pdat, size = 0.3,
                  legend.labs = unique(pdat$CDC3_group),
                  surv.median.line = "none",
                  pval = TRUE,
                  conf.int = FALSE,
                  palette = c("grey", "#00726A"),
                  title = paste(title_prefix, "Overall Survival"))
  print(p)
  
  # ---- 统计检验 ----
  logrank_test <- survdiff(Surv(OS.time, BR2) ~ CDC3_group, data = pdat)
  logrank_pval <- signif(1 - pchisq(logrank_test$chisq, df = 1), digits = 4)
  
  coxregre <- coxph(Surv(OS.time, BR2) ~ CDC3_group, data = pdat)
  hr <- signif(summary(coxregre)$coef[1, "exp(coef)"], digits = 2)
  ci.up <- signif(summary(coxregre)$conf.int[1, 4], digits = 2)
  ci.down <- signif(summary(coxregre)$conf.int[1, 3], digits = 2)
  
  # ---- 输出结果 ----
  tab <- table(pdat$CDC3_group)
  cat("\n===", title_prefix, "Survival Analysis Results ===\n")
  cat("Sample size:", nrow(pdat), "\n")
  cat("High group n:", tab["high"], "| Low group n:", tab["low"], "\n")
  cat("Hazard Ratio (HR):", hr, "\n")
  cat("95% CI:", ci.down, "-", ci.up, "\n")
  cat("Log-rank P-value:", logrank_pval, "\n")
  cat("Cox PH P-value:", signif(summary(coxregre)$coefficients[1, 5], digits = 4), "\n")
  
  return(list(plot = p, fit = fit, cox = coxregre, data = pdat))
}

# ============ 3. 分别运行 ============
pre_results  <- run_survival_analysis(gexp_pre,  pheno_pre,  "Pre-treatment")
tiff("~/202311_tnbc_caf_05_20240521/GDC.download/KM_PRJEB.tiff",
     width = 105, height = 110, units = "mm", res = 600,   # 单栏或1.5栏尺寸
     compression = "lzw")
print(pre_results$plot)
dev.off()

post_results <- run_survival_analysis(gexp_post, pheno_post, "Post-treatment")
pre_results$plot
#ROC###
pre_pdat <- filter(pdat,pdat$Treatment=="PRE")
library(pROC)
pre_pdat=pre_pdat %>% mutate(BR2 = ifelse(response_NR =="R", 1, 0)) 

proc <- roc(pre_pdat$BR2, pre_pdat$IFNG_signature,
            smoothed = T,
            ci = T, ci.alpha = 0.9, stratified = F,
            plot = T, auc.polygon = T, max.auc.polygon = T, grid = T, print.auc = T, show.thres = T)

proc_auc <- round(as.numeric(proc$ci),2)

proc_ncases = length(proc$cases) + length(proc$controls)

# function to extract data for plotting ROC curves
roc2df <- function(p) {
  
  tibble(TPR = p$sensitivities,
         FPR = 1 - p$specificities)
  
}
roc2df(proc) %>%
  distinct %>%
  dplyr::arrange(TPR) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_line(size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 1, size = 1, colour = "grey50") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "Melanoma-PRJEB23709") +
  annotate(geom = "text", x = 1, y = 0.05, label = paste("n = ", proc_ncases, "\nAUC ", proc_auc[2], collapse = "", sep = ""), size = 5, hjust = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5))

verification::roc.area(as.numeric(as.vector(pre_pdat$BR2)),proc$predictor)

