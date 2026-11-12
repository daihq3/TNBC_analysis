library(readr)
library(dplyr)
library(ggplot2)

# ── 1. Read MISTY importance output ────────────────────────────────────────────
imp_file <- "~/residual spatial/mistyr/MISTY_RESULT/cm_CID4465_c2l/importances_CD8.CXCL13_para_15.txt"
imp_all <- readr::read_csv(imp_file, show_col_types = FALSE)

# Standardize column names
if ("target" %in% colnames(imp_all)) {
  imp_all <- imp_all %>%
    dplyr::rename(predictor = target)
}

# ── 2. Define myeloid predictors of interest ──────────────────────────────────
myeloid_predictors <- c(
  "Macro.C1QC", "Macro.CXCL10", "Macro.FOLR2", "Macro.SPP1",
  "Mono.FCN1", "cDC1.CLEC9A", "cDC2.CD1C", "cDC3.LAMP3", "pDC.LILRA4"
)

# Clean predictor labels: Macro.C1QC → Mac-C1QC, cDC3.LAMP3 → cDC3-LAMP3, etc.
clean_label <- function(x) {
  x <- gsub("^Macro\\.", "Mac-", x)
  x <- gsub("^Mono\\.", "Mono-", x)
  x <- gsub("^cDC", "cDC", x)
  x <- gsub("^pDC\\.", "pDC-", x)
  x <- gsub("\\.", "-", x)
  return(x)
}

imp_myeloid <- imp_all %>%
  dplyr::filter(predictor %in% myeloid_predictors) %>%
  dplyr::mutate(
    predictor_clean = clean_label(predictor),
    # ensure consistent ordering: keep original factor order after sorting by importance
    predictor_clean = factor(predictor_clean, levels = predictor_clean[order(imp)])
  )

# ── 3. Identify top predictor ─────────────────────────────────────────────────
top_predictor <- imp_myeloid %>%
  dplyr::slice_max(order_by = imp, n = 1) %>%
  dplyr::pull(predictor_clean)

message("Top myeloid predictor: ", top_predictor)

# ── 4. Single-row heatmap ─────────────────────────────────────────────────────
p <- ggplot(imp_myeloid,
            aes(x = predictor_clean, y = "CD8-CXCL13", fill = imp)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # highlight the top predictor with a black border
  geom_tile(
    data = imp_myeloid %>% dplyr::filter(predictor_clean == top_predictor),
    aes(x = predictor_clean, y = "CD8-CXCL13"),
    color = "black", linewidth = 1.5, fill = NA
  ) +
  scale_fill_gradient2(
    low  = "#2166AC",   # blue
    mid  = "#F7F7F7",   # white
    high = "#B2182B",   # red
    midpoint = median(imp_myeloid$imp),
    name = "Median\nimportance",
    limits = c(0, max(imp_myeloid$imp)),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "MISTY spatial dependency: CD8-CXCL13",
    subtitle = paste0("Myeloid predictors (CID4465); top predictor: ", top_predictor),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, colour = "grey40", hjust = 0.5),
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y     = element_text(size = 12, face = "bold"),
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    panel.grid      = element_blank(),
    plot.margin     = margin(10, 20, 10, 10)
  )

print(p)

# ── 5. Save ────────────────────────────────────────────────────────────────────
out_dir <- "D:/202509/FOXP1/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "MISTY_CD8CXCL13_myeloid_predictors_heatmap.pdf"),
       p, width = 8, height = 3.5, useDingbats = FALSE)
ggsave(file.path(out_dir, "MISTY_CD8CXCL13_myeloid_predictors_heatmap.png"),
       p, width = 8, height = 3.5, dpi = 300)

# ── 6. Output sorted table ─────────────────────────────────────────────────────
out_df <- imp_myeloid %>%
  dplyr::arrange(desc(imp)) %>%
  dplyr::select(predictor, predictor_clean, imp)
print(out_df)

write.csv(out_df,
          file.path(out_dir, "MISTY_CD8CXCL13_myeloid_predictors.csv"),
          row.names = FALSE)

message("Done. Files saved to: ", out_dir)
