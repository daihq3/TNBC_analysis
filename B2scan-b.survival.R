#load(https://github.com/tkarn/G9-HTG/tree/main/input/valid_datasets_HTG.RData)

#Molecular adaptation to neoadjuvant immunotherapy in triple-negative breast cancer
load("~/valid_datasets_HTG.RData")
markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")
CD8.CXCL13 <- markers_fine_celltype %>%
  filter(cluster == "CD8-CXCL13" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)
geneset_list=list(IFNG_signature=CD8.CXCL13$gene)
YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")
# saveRDS(YS.combined_markers,"YS.combined_markers.rds")
cdc.lamp3 <- YS.combined_markers[YS.combined_markers$cluster=="cDC3-LAMP3",]
selected_genes <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)
selected_genes
geneset_list=list(IFNG_signature=selected_genes$gene)

gexp.df <-  scanb.htg.Rseq

# gexp.df <- t(gexp.df)
# colnames(gexp.df) <- gexp.df[1,]
gexp.df <- as.matrix(gexp.df)
# gexp.df <- gexp.df[-1,]




selected_phenos <- scanb.htg.cli.dat
pdat <-
  selected_phenos %>%
  as.data.frame() 
pdat$BR2=pdat$OS_event
pdat$OS.time <- pdat$OS_months/12
choose_samples1 <- pdat %>%
    filter(PAM50 == "Basal") %>%
    pull(X_Sample_title)

expr_TNBC <- gexp.df[, choose_samples1]
expr_mat <- as.matrix(expr_TNBC)
expr_mat[is.na(expr_mat)] <- 0
storage.mode(expr_mat) <- "numeric"
gsvapar <- gsvaParam(expr_mat, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)
pdat.basal <- pdat[pdat$X_Sample_title %in%choose_samples1,]
pdat.basal$inf_sig_score <- t(gsva_es)


pdat.basal$inf_sig_score <- as.numeric(pdat.basal$inf_sig_score)

pdat.basal$OS.time <- pdat.basal$OS_years


res.cut <- surv_cutpoint(pdat.basal, time = "OS.time", 
                         event = "BR2", 
                         variables = "inf_sig_score")
res.cat <- surv_categorize(res.cut)
pdat.basal$CDC3_group<- res.cat[,"inf_sig_score"] 
pdat.basal$CDC3_group = factor(pdat.basal$CDC3_group, levels = c("low","high")) 
fit = survfit(Surv(OS.time,BR2) ~ CDC3_group, data=pdat.basal)


p = ggsurvplot(fit, pdat.basal, size=0.3, vlegend.labs=unique(pdat.basal$CDC3_group),
               surv.median.line="none", pval=T, conf.int=F,
               palette=c("grey","#00726A"))
### (1) 各组人数、中位生存时间、log-rank P —— surv_summary 一行搞定 ###
surv_summary(fit, data = pdat.basal)
# 输出表格中包含：n（各组起始人数）、median（中位生存时间）等

### (2) HR 和 95%CI —— Cox回归 ###
cox_res <- coxph(Surv(OS.time,BR2) ~ CDC3_group, data = pdat.basal)
summary(cox_res)
# 其中 coef = log(HR)，HR = exp(coef)
# high vs low 的 HR 以 low 组为参照（因为factor levels设了low在前）

# 只提取干净的HR表格
hr_table <- summary(cox_res)$conf.int   # HR、95%CI上下限
hr_table

### (3) log-rank P 值（显式计算，可写入论文）###
sd <- survdiff(Surv(OS.time,BR2) ~ CDC3_group, data=pdat.basal)
p_value <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_value


# 组装标注文本
hr <- round(hr_table[1, 1], 2)
ci_low <- round(hr_table[1, 3], 2)
ci_high <- round(hr_table[1, 4], 2)
n_low <- table(pdat.basal$CDC3_group)["low"]
n_high <- table(pdat.basal$CDC3_group)["high"]

label_text <- paste0("HR = ", hr, " (", ci_low, "-", ci_high, ")\n",
                     "n(low) = ", n_low, ", n(high) = ", n_high, "\n",
                     "Log-rank P = ", signif(p_value, 3))
label_text
p$plot + 
  annotate("text", x = max(pdat.basal$OS.time) * 0.6, 
           y = 0.9, label = label_text, hjust = 0, size = 3.5)

tiff("~/202311_tnbc_caf_05_20240521/GDC.download/KM_TNBC.tiff",
     width = 120, height = 110, units = "mm", res = 600,   # 单栏或1.5栏尺寸
     compression = "lzw")
print(p)
dev.off()



