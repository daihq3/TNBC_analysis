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

library(tidyr)
pd <- pd %>% inner_join(OS,by = "PATIENT_ID")  
pd <- pd %>%  
  mutate(`OS.time`=OS_MONTHS/12) %>%  
  separate(OS_STATUS,into = c("OS","STATUS"),sep = ":")  
pd$Type <- factor(pd$Type,levels=c("TNBC","HER2","Luminal A","Luminal B"))
df <- as.data.frame(t(exprSet)) %>%  
  rownames_to_column("PATIENT_ID") %>% inner_join(pd,by = "PATIENT_ID")  
head(df[,1:6])  


markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")
selected_genes <- markers_fine_celltype %>%
    filter(cluster == "CD8-CXCL13" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
    dplyr::select(gene)

library(dplyr)
geneset_list=list(IFNG_signature=selected_genes$gene)

choose_samples1 <- pd %>%
  filter(Type == "TNBC") %>%
  pull(PATIENT_ID)

expr_TNBC <- expr[,choose_samples1]

expr_mat <- as.matrix(expr_TNBC)
expr_mat[is.na(expr_mat)] <- 0 #for HER2
storage.mode(expr_mat) <- "numeric"
gsvapar <- gsvaParam(expr_mat, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)
library(data.table)

Metabric.tnbc.meta <- fread("~/202311_tnbc_caf_05_20240521/GDC.download/brca_metabric/data_clinical_patient.txt")
Metabric.tnbc.meta <- Metabric.tnbc.meta %>% dplyr::filter(PATIENT_ID %in% choose_samples1)
dim(Metabric.tnbc.meta) #[1] 334  24
dim(Metabric.tnbc.meta)

Metabric.tnbc.meta$inf_sig_score <- t(gsva_es)
Metabric.tnbc.meta$inf_sig_score <- as.numeric(Metabric.tnbc.meta$inf_sig_score)
Metabric.tnbc.meta$OS_STATUS <- gsub(":.*", "", Metabric.tnbc.meta$OS_STATUS)
Metabric.tnbc.meta$OS_STATUS <- as.numeric(Metabric.tnbc.meta$OS_STATUS)
Metabric.tnbc.meta$OS_MONTHS <- as.numeric(Metabric.tnbc.meta$OS_MONTHS)
Metabric.tnbc.meta$OS.time <- Metabric.tnbc.meta$OS_MONTHS/12

res.cut <- surv_cutpoint(Metabric.tnbc.meta, time = "OS.time", 
                         event = "OS_STATUS", 
                         variables = "inf_sig_score")
res.cat <- surv_categorize(res.cut)
Metabric.tnbc.meta$CD8CXCL13_group<- res.cat[,"inf_sig_score"] 
Metabric.tnbc.meta$CD8CXCL13_group = factor(Metabric.tnbc.meta$CD8CXCL13_group, levels = c("low","high"))

fit = survfit(Surv(OS.time,OS_STATUS) ~ CD8CXCL13_group, data=Metabric.tnbc.meta)
p = ggsurvplot(fit, Metabric.tnbc.meta, size=0.3, vlegend.labs=unique(Metabric.tnbc.meta$CD8CXCL13_group),
               surv.median.line="none", pval=T, conf.int=F,
               palette=c("grey","#00726A"))


### (1) 各组人数、中位生存时间、log-rank P —— surv_summary 一行搞定 ###
surv_summary(fit, data = Metabric.tnbc.meta)
# 输出表格中包含：n（各组起始人数）、median（中位生存时间）等

### (2) HR 和 95%CI —— Cox回归 ###
cox_res <- coxph(Surv(OS.time, OS_STATUS) ~ CD8CXCL13_group, data = Metabric.tnbc.meta)
summary(cox_res)
# 其中 coef = log(HR)，HR = exp(coef)
# high vs low 的 HR 以 low 组为参照（因为factor levels设了low在前）

# 只提取干净的HR表格
hr_table <- summary(cox_res)$conf.int   # HR、95%CI上下限
hr_table

### (3) log-rank P 值（显式计算，可写入论文）###
sd <- survdiff(Surv(OS.time, OS_STATUS) ~ CD8CXCL13_group, data = Metabric.tnbc.meta)
p_value <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_value


# 组装标注文本
hr <- round(hr_table[1, 1], 2)
ci_low <- round(hr_table[1, 3], 2)
ci_high <- round(hr_table[1, 4], 2)
n_low <- table(Metabric.tnbc.meta$CD8CXCL13_group)["low"]
n_high <- table(Metabric.tnbc.meta$CD8CXCL13_group)["high"]

label_text <- paste0("HR = ", hr, " (", ci_low, "-", ci_high, ")\n",
                     "n(low) = ", n_low, ", n(high) = ", n_high, "\n",
                     "Log-rank P = ", signif(p_value, 3))
label_text
p$plot + 
  annotate("text", x = max(Metabric.tnbc.meta$OS.time) * 0.6, 
           y = 0.9, label = label_text, hjust = 0, size = 3.5)

tiff("~/202311_tnbc_caf_05_20240521/GDC.download/KM_TNBC.tiff",
     width = 120, height = 110, units = "mm", res = 600,   # 单栏或1.5栏尺寸
     compression = "lzw")
print(p)
dev.off()

YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")


selected_genes <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)
#120
library(dplyr)
geneset_list=list(IFNG_signature=selected_genes$gene)

-----------------------------no-treatment with chemo-------------------------------------


-----------------------------forest plot---
-------------------------------------

library(dplyr)

pd_filtered <- Metabric.tnbc.meta %>%
    left_join(pd, by = "PATIENT_ID")

pdat <-pd_filtered

res.cut <- surv_cutpoint(pdat, time = "OS.time.x", 
                         event = "OS_STATUS", 
                         variables = "inf_sig_score")
res.cat <- surv_categorize(res.cut)
pdat$CDC3 <- res.cat[,"inf_sig_score"] 
pdat$CDC3 = factor(pdat$CDC3)
pdat$CDC3 = factor(pdat$CDC3, levels = c("low","high"))


clin_merge1 <- pdat %>%
    dplyr::filter(
        !is.na(AGE_AT_DIAGNOSIS.x),
        !is.na(LYMPH_NODES_EXAMINED_POSITIVE.x),
        !is.na(GRADE),
        !is.na(CHEMOTHERAPY.x),
        !is.na(TUMOR_SIZE),
        !is.na(OS.time.x),
        !is.na(OS_STATUS),
        !is.na(CDC3),
        HISTOLOGICAL_SUBTYPE.x != "" | is.na(HISTOLOGICAL_SUBTYPE.x)  # 空字符串视情况
    ) %>%
    dplyr::mutate(
        AGE_AT_DIAGNOSIS = as.numeric(AGE_AT_DIAGNOSIS.x),
        Age.group = ifelse(AGE_AT_DIAGNOSIS >= 50, ">=50", "<50"),
        
        LYMPH_NODES_EXAMINED_POSITIVE = as.numeric(LYMPH_NODES_EXAMINED_POSITIVE.x),
        
        # Grade 分组：以 grade 1 为参照
        Grade = as.numeric(GRADE),
        Grade_group = factor(ifelse(GRADE %in% c(1, 2), "Grade 1/2", "Grade 3"),
                             levels = c("Grade 3", "Grade 1/2")),  # reference = Grade 3
        
        nodal.status = factor(ifelse(LYMPH_NODES_EXAMINED_POSITIVE == 0, "No", "Yes")),
        
        # 化疗：YES vs NO（reference = NO）
        Chemo_group = factor(
            ifelse(is.na(CHEMOTHERAPY.x), NA,
                   ifelse(toupper(trimws(CHEMOTHERAPY.x)) == "YES", "Yes", "No")),
            levels = c("No", "Yes")   # reference = No
        ),
        
        # 组织学分组
        Histo_group = ifelse(HISTOLOGICAL_SUBTYPE.x == "Ductal/NST", "Ductal/NST", "Other"),
        Histo_group = factor(Histo_group, levels = c("Ductal/NST", "Other")),
        
        # 肿瘤大小分组
        tumor_size_group = factor(
            case_when(
                TUMOR_SIZE < 20  ~ "Small (<2cm)",
                TUMOR_SIZE >= 20 ~ "Large (>=2cm)"),
            levels = c("Small (<2cm)", "Large (>=2cm)")
        )
    )
coxregre_multivariable <- coxph(
    Surv(OS.time.x, OS_STATUS) ~ CDC3 + Age.group + nodal.status + tumor_size_group + Histo_group + Grade_group + Chemo_group,
    data = clin_merge1
)
coxregre_multivariable
summary(coxregre_multivariable)

#按照数值改
library(ggplot2)
library(dplyr)

# 按 coxregre_multivariable 的结果更新数据
forest_data <- data.frame(
  variable = c("CDC3 (high)", "Age >= 50", "Nodal Status (Yes)",
               "Tumor Size (>=2cm)", "Histo_group (Other)",
               "Grade 1/2", "Chemo (Yes)"),
  # exp(coef)
  HR = c(0.56058, 1.52830, 1.88249, 1.47643, 0.93540, 0.99954, 0.73700),
  # 95% CI
  lower = c(0.4059, 1.0619, 1.1618, 1.0145, 0.5917, 0.6175, 0.4491),
  upper = c(0.7743, 2.1996, 3.0502, 2.1488, 1.4787, 1.6180, 1.2095),
  # Pr(>|z|)
  p_value = c(0.000443, 0.022415, 0.010199, 0.041851, 0.775037, 0.998507, 0.227308)
)

forest_data <- forest_data %>%
  mutate(variable = factor(variable, levels = variable[order(HR)]))

library(ggplot2)
library(dplyr)

forest_data <- data.frame(
  variable = c("CDC3 (high)", "Age >= 50", "Nodal Status (Yes)",
               "Tumor Size (>=2cm)", "Histo_group (Other)",
               "Grade 1/2", "Chemo (Yes)"),
  HR = c(0.56058, 1.52830, 1.88249, 1.47643, 0.93540, 0.99954, 0.73700),
  lower = c(0.4059, 1.0619, 1.1618, 1.0145, 0.5917, 0.6175, 0.4491),
  upper = c(0.7743, 2.1996, 3.0502, 2.1488, 1.4787, 1.6180, 1.2095),
  p_value = c(0.000443, 0.022415, 0.010199, 0.041851, 0.775037, 0.998507, 0.227308)
)

forest_data <- forest_data %>%
  mutate(variable = factor(variable, levels = variable[order(HR)]))

p <- ggplot(forest_data, aes(x = HR, y = variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "grey50") +
  geom_point(aes(color = p_value < 0.05), size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  # 直接显示原始 p 值，保留足够位数
  geom_text(aes(label = ifelse(p_value < 0.05,
                               paste0("p=", format(p_value, scientific = FALSE,
                                                   digits = 3)), "")),
            vjust = -1.5, size = 3) +
  scale_color_manual(values = c("TRUE" = "#D32F2F", "FALSE" = "black")) +
  labs(title = "Multivariate Cox Regression for Overall Survival",
       x = "Hazard Ratio", y = "") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())

print(p)



-------
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

library(tidyr)

pd <- pd %>% inner_join(OS,by = "PATIENT_ID")  
pd <- pd %>%  
    mutate(`OS.time`=OS_MONTHS/12) %>%  
    separate(OS_STATUS,into = c("OS","STATUS"),sep = ":")  
pd$Type <- factor(pd$Type,levels=c("TNBC","HER2","Luminal A","Luminal B"))
df <- as.data.frame(t(exprSet)) %>%  
    rownames_to_column("PATIENT_ID") %>% inner_join(pd,by = "PATIENT_ID")  
head(df[,1:6])  

choose_samples1 <- pd %>%
    filter(Type == "TNBC") %>%
    pull(PATIENT_ID)

YS.combined_markers <- readRDS("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/LAMP3-CD8/YS.combined_markers.rds")


selected_genes <- YS.combined_markers %>%
  filter(cluster == "cDC3-LAMP3" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene)
#120
library(dplyr)
geneset_list=list(IFNG_signature=selected_genes$gene)

expr_TNBC <- expr[,choose_samples1]

expr_mat <- as.matrix(expr_TNBC)
expr_mat[is.na(expr_mat)] <- 0 #for HER2
storage.mode(expr_mat) <- "numeric"
gsvapar <- gsvaParam(expr_mat, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)
library(data.table)

Metabric.tnbc.meta <- fread("~/202311_tnbc_caf_05_20240521/GDC.download/brca_metabric/data_clinical_patient.txt")
Metabric.tnbc.meta <- Metabric.tnbc.meta %>% dplyr::filter(PATIENT_ID %in% choose_samples1)
dim(Metabric.tnbc.meta) #[1] 334  24
dim(Metabric.tnbc.meta)

Metabric.tnbc.meta$inf_sig_score <- t(gsva_es)


markers_fine_celltype <- read_csv("~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/markers_fine_celltype.csv")
selected_genes <- markers_fine_celltype %>%
    filter(cluster == "CD8-CXCL13" & avg_log2FC >= 1 & p_val_adj < 0.05) %>%
    dplyr::select(gene)

library(dplyr)
geneset_list=list(IFNG_signature=selected_genes$gene)

choose_samples1 <- pd %>%
    filter(Type == "TNBC") %>%
    pull(PATIENT_ID)

expr_TNBC <- expr[,choose_samples1]

expr_mat <- as.matrix(expr_TNBC)
expr_mat[is.na(expr_mat)] <- 0 #for HER2
storage.mode(expr_mat) <- "numeric"
gsvapar <- gsvaParam(expr_mat, geneset_list, maxDiff = TRUE)
gsva_es <- gsva(gsvapar)
library(data.table)

Metabric.tnbc.meta$inf_sig_score2 <- t(gsva_es)




pd_filtered <- Metabric.tnbc.meta %>%
    left_join(pd, by = "PATIENT_ID")

pdat <-pd_filtered


pdat$OS_STATUS <- gsub(":.*", "", pdat$OS_STATUS)
pdat$OS_STATUS <- as.numeric(pdat$OS_STATUS)
pdat$OS_MONTHS <- as.numeric(pdat$OS_MONTHS)
pdat$OS.time <- pdat$OS_MONTHS.y/12


# ---- CDC3 signature: median split ----
pdat$CDC3 <- ifelse(pdat$inf_sig_score > median(pdat$inf_sig_score, na.rm = TRUE),
                    "high", "low")
pdat$CDC3 <- factor(pdat$CDC3, levels = c("low", "high"))

# ---- CD8-CXCL13 signature: median split ----
pdat$CD8CXCL13 <- ifelse(pdat$inf_sig_score2 > median(pdat$inf_sig_score2, na.rm = TRUE),
                         "high", "low")
pdat$CD8CXCL13 <- factor(pdat$CD8CXCL13, levels = c("low", "high"))
# 检查各组的样本数，应该大致各占一半
table(pdat$CDC3)
table(pdat$CD8CXCL13)
clin_merge2 <- pdat %>%
    dplyr::mutate(
        CD8_CDC3_Group = case_when(
            CD8CXCL13 == "high" & CDC3 == "high" ~ "CD8_High/CDC3_High",
            CD8CXCL13 == "high" & CDC3 == "low"  ~ "CD8_High/CDC3_Low",
            CD8CXCL13 == "low"  & CDC3 == "high" ~ "CD8_Low/CDC3_High",
            CD8CXCL13 == "low"  & CDC3 == "low"  ~ "CD8_Low/CDC3_Low"
        )
    ) %>%
    # 设置因子的参考水平
    # 通常将 "Double Low" (双低) 设为 Reference (参考组)
    mutate(CD8_CDC3_Group = factor(CD8_CDC3_Group, 
                                   levels = c("CD8_Low/CDC3_Low",
                                              "CD8_Low/CDC3_High",
                                              "CD8_High/CDC3_Low",
                                              "CD8_High/CDC3_High")))
fit = survfit(Surv(OS.time,OS_STATUS) ~ CD8_CDC3_Group, data=clin_merge2)
p = ggsurvplot(fit, clin_merge2, size=0.3, vlegend.labs=unique(pdat$CD8_CDC3_Group),
                surv.median.line="none", pval=T, conf.int=F,
                palette=c("#e31a1b","#cbcbcb","#7f7f7f","#367db7"))
 
p

### (1) 各组人数、中位生存时间、log-rank P —— surv_summary 一行搞定 ###
surv_summary(fit, data = clin_merge2)
# 输出表格中包含：n（各组起始人数）、median（中位生存时间）等

### (2) HR 和 95%CI —— Cox回归 ###
cox_res <- coxph(Surv(OS.time, OS_STATUS) ~ CD8_CDC3_Group, data = clin_merge2)
summary(cox_res)
# 其中 coef = log(HR)，HR = exp(coef)
# high vs low 的 HR 以 low 组为参照（因为factor levels设了low在前）

# 只提取干净的HR表格
hr_table <- summary(cox_res)$conf.int   # HR、95%CI上下限
hr_table

### (3) log-rank P 值（显式计算，可写入论文）###
sd <- survdiff(Surv(OS.time, OS_STATUS) ~ CD8_CDC3_Group, data = clin_merge2)
p_value <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_value

clin_merge1 <- pdat %>%
    dplyr::filter(
        !is.na(AGE_AT_DIAGNOSIS.x),
        !is.na(LYMPH_NODES_EXAMINED_POSITIVE.x),
        !is.na(GRADE),
        !is.na(CHEMOTHERAPY.x),
        !is.na(TUMOR_SIZE),
        !is.na(OS.time),
        !is.na(OS_STATUS),
        !is.na(CDC3),
        HISTOLOGICAL_SUBTYPE.x != "" | is.na(HISTOLOGICAL_SUBTYPE.x)
    ) %>%
    dplyr::mutate(
        AGE_AT_DIAGNOSIS = as.numeric(AGE_AT_DIAGNOSIS.x),
        Age.group = ifelse(AGE_AT_DIAGNOSIS >= 50, ">=50", "<50"),
        
        LYMPH_NODES_EXAMINED_POSITIVE = as.numeric(LYMPH_NODES_EXAMINED_POSITIVE.x),
        
        Grade = as.numeric(GRADE),
        Grade_group = factor(ifelse(GRADE %in% c(1, 2), "Grade 1/2", "Grade 3"),
                             levels = c("Grade 3", "Grade 1/2")),
        
        nodal.status = factor(ifelse(LYMPH_NODES_EXAMINED_POSITIVE == 0, "No", "Yes")),
        
        Chemo_group = factor(
            ifelse(is.na(CHEMOTHERAPY.x), NA,
                   ifelse(toupper(trimws(CHEMOTHERAPY.x)) == "YES", "Yes", "No")),
            levels = c("No", "Yes")
        ),
        
        Histo_group = factor(
            ifelse(HISTOLOGICAL_SUBTYPE.x == "Ductal/NST", "Ductal/NST", "Other"),
            levels = c("Ductal/NST", "Other")
        ),
        
        tumor_size_group = factor(
            case_when(
                TUMOR_SIZE < 20  ~ "Small (<2cm)",
                TUMOR_SIZE >= 20 ~ "Large (>=2cm)"
            ),
            levels = c("Small (<2cm)", "Large (>=2cm)")
        ),   # ← 原代码这里缺了两个右括号和逗号
        
        CD8_CDC3_Group = case_when(
            CD8CXCL13 == "high" & CDC3 == "high" ~ "CD8_High/CDC3_High",
            CD8CXCL13 == "high" & CDC3 == "low"  ~ "CD8_High/CDC3_Low",
            CD8CXCL13 == "low"  & CDC3 == "high" ~ "CD8_Low/CDC3_High",
            CD8CXCL13 == "low"  & CDC3 == "low"  ~ "CD8_Low/CDC3_Low"
        )
    ) %>%
    mutate(CD8_CDC3_Group = factor(CD8_CDC3_Group,
                                   levels = c("CD8_Low/CDC3_Low",
                                              "CD8_Low/CDC3_High",
                                              "CD8_High/CDC3_Low",
                                              "CD8_High/CDC3_High")))

coxregre_multivariable <- coxph(
    Surv(OS.time.x, OS_STATUS) ~ CD8_CDC3_Group + Age.group + nodal.status + tumor_size_group + Histo_group + Grade_group + Chemo_group,
    data = clin_merge1
)
coxregre_multivariable
summary(coxregre_multivariable)


--------去掉chemo-both---
library(ggplot2)
library(dplyr)

forest_data <- data.frame(
  variable = c("CD8_High/CDC3_High", "CD8_Low/CDC3_High", "CD8_High/CDC3_Low",
               "Age >= 50", "Nodal Status (Yes)",
               "Tumor Size (>=2cm)", "Histo_group (Other)",
               "Grade 1/2"),
  HR = c(0.507581, 0.577382, 0.652597, 1.622905, 1.547254, 
         1.471660, 1.006814, 1.070840),
  lower = c(0.3477, 0.3336, 0.3873, 1.1258, 1.1102, 
            1.0111, 0.6325, 0.6695),
  upper = c(0.7409, 0.9992, 1.0996, 2.3394, 2.1564, 
            2.1419, 1.6027, 1.7129),
  p_value = c(0.000442, 0.049666, 0.108860, 0.009451, 0.009963, 
              0.043611, 0.977160, 0.775196)
)

forest_data <- forest_data %>%
  mutate(variable = factor(variable, levels = variable[order(HR)]))

p <- ggplot(forest_data, aes(x = HR, y = variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "grey50") +
  geom_point(aes(color = p_value < 0.05), size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  geom_text(aes(label = ifelse(p_value < 0.05,
                               paste0("p=", format(p_value, scientific = FALSE,
                                                   digits = 3)), "")),
            vjust = -1.5, size = 3) +
  scale_color_manual(values = c("TRUE" = "#D32F2F", "FALSE" = "black")) +
  labs(title = "Multivariate Cox Regression for Overall Survival",
       x = "Hazard Ratio", y = "") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())

print(p)

ggsave("forest_plot_multivariable.pdf", p, width = 7, height = 5)






