library(Seurat)
library(Nebulosa)
library(ggplot2)
library(dplyr)

# ── Setup ─────────────────────────────────────────────────────────────────────
# Tcell     : Combined_TNBC_Lambrecht T cell subset (has UMAP)
# Tcell.filt: FOXP1 project T cell subset (has CD8-CXCL13 in cellstate)

# ── 1. Get barcodes of CD8-CXCL13 cells from Tcell.filt ─────────────────────
cd8cx_barcodes <- rownames(Tcell.filt@meta.data[Tcell.filt$cellstate == "CD8-CXCL13", ])
message("CD8-CXCL13 cells in Tcell.filt: ", length(cd8cx_barcodes))

# ── 2. Build per-cluster CD8+ signatures from fine-cell-type markers ─────────
# For each CD8+ cluster, take top 100 marker genes ranked by avg_log2FC,
# then compute a separate module score on Tcell.Lam.

cd8_clusters <- c("CD8-GZMK", "CD8-CCL5", "CD8-CXCL13", "CD8-gdT-TRDC")

# Build a named list of signatures
cluster_sigs <- lapply(cd8_clusters, function(cl) {
  markers_fine_celltype %>%
    dplyr::filter(cluster == cl) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::pull(gene) %>%
    unique()
})
names(cluster_sigs) <- cd8_clusters

# Print signature sizes
for (cl in cd8_clusters) {
  message(cl, " signature genes: ", length(cluster_sigs[[cl]]))
  message("  Top 10: ", paste(head(cluster_sigs[[cl]], 10), collapse = ", "))
}

# ── 3. Score all cells in Tcell.Lam with each per-cluster signature ───────────
for (cl in cd8_clusters) {
  score_name <- make.names(paste0("score_", cl))
  Tcell.Lam <- AddModuleScore(
    Tcell.Lam,
    features = list(cluster_sigs[[cl]]),
    name     = score_name,
    assay    = "RNA",
    seed     = 42
  )
  # AddModuleScore appends "1"; create positive-clipped version
  score_col <- paste0(score_name, "1")
  pos_col   <- paste0(score_name, "_pos")
  Tcell.Lam[[pos_col]] <- pmax(Tcell.Lam[[score_col]][, 1], 0)

  message(cl, " score range: ",
          round(min(Tcell.Lam[[score_col]][, 1]), 3), " to ",
          round(max(Tcell.Lam[[score_col]][, 1]), 3))
}

# Use harmonized object name for downstream plotting
Tcell <- Tcell.Lam

# ── Option A: Nebulosa density for each per-cluster signature ─────────────────
density_plots <- list()

for (cl in cd8_clusters) {
  score_name <- make.names(paste0("score_", cl))
  pos_col    <- paste0(score_name, "_pos")

  p <- plot_density(
    Tcell,
    features  = pos_col,
    reduction = "umap",
    method    = "ks"
  ) +
    labs(title = paste0(cl, " signature density on Tcell.Lam UMAP")) +
    theme_classic(base_size = 11) +
    theme(
      plot.title  = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8)
    )

  density_plots[[cl]] <- p
  print(p)
}

# Combine all density plots into one figure
library(patchwork)
p_all_density <- wrap_plots(density_plots, ncol = 2) +
  plot_annotation(title = "Per-cluster CD8+ signature densities on Tcell.Lam UMAP")
print(p_all_density)

# Message positive-score cell counts
for (cl in cd8_clusters) {
  score_name <- make.names(paste0("score_", cl))
  pos_col    <- paste0(score_name, "_pos")
  message("Cells with positive ", cl, " score: ", sum(Tcell[[pos_col]][, 1] > 0))
}

# ── Option B: UMAP with CD8-CXCL13 highlighted (for comparison) ──────────────
umap_df <- as.data.frame(Tcell@reductions$umap@cell.embeddings)
# Normalize column names (handle UMAP_1/UMAP_2 or umap_1/umap_2)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$highlight <- ifelse(colnames(Tcell) %in% cd8cx_barcodes,
                             "CD8-CXCL13", "Other")
umap_df$highlight <- factor(umap_df$highlight, levels = c("Other", "CD8-CXCL13"))

p_highlight <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = highlight)) +
  geom_point(data = subset(umap_df, highlight == "Other"),
             aes(x = UMAP_1, y = UMAP_2),
             size = 0.3, alpha = 0.3, color = "lightgrey") +
  geom_point(data = subset(umap_df, highlight == "CD8-CXCL13"),
             aes(x = UMAP_1, y = UMAP_2),
             size = 0.5, alpha = 0.8, color = "#1a7a1a") +
  labs(
    title = "CD8-CXCL13 (Tcell.filt) on Tcell UMAP",
    x = "UMAP 1", y = "UMAP 2",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(p_highlight)

# ── Save ──────────────────────────────────────────────────────────────────────
out_dir <- "D:/202509/FOXP1/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save combined density figure
ggsave(file.path(out_dir, "umap_cd8_all_clusters_density.pdf"),
       p_all_density, width = 10, height = 10, useDingbats = FALSE)

# Save individual cluster density plots
for (cl in cd8_clusters) {
  safe_name <- make.names(cl)
  ggsave(file.path(out_dir, paste0("umap_", safe_name, "_density.pdf")),
         density_plots[[cl]], width = 6, height = 5, useDingbats = FALSE)
}

ggsave(file.path(out_dir, "umap_cd8cxcl13_highlight.pdf"), p_highlight,
       width = 6, height = 5, useDingbats = FALSE)
message("Saved per-cluster density and highlight UMAPs.")
