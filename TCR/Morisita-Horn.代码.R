
# ── 13. Clone overlap between cell types (Morisita-Horn index) ────────────────
# Build cell type × clonotype count matrix
ct_clone_mat <- seu@meta.data %>%
  dplyr::filter(both_chains == "both") %>%
  dplyr::count(cell_type_fine, clonotype_id) %>%
  tidyr::pivot_wider(names_from = cell_type_fine,
                     values_from = n,
                     values_fill = 0) %>%
  tibble::column_to_rownames("clonotype_id") %>%
  as.matrix()

# Morisita-Horn similarity between cell types
morisita_horn <- function(x, y) {
  xy_sum  <- sum(x * y)
  x2_sum  <- sum(x^2)
  y2_sum  <- sum(y^2)
  2 * xy_sum / ((x2_sum / sum(x)^2 + y2_sum / sum(y)^2) * sum(x) * sum(y))
}

ct_names <- colnames(ct_clone_mat)
n_ct     <- length(ct_names)
overlap_mat <- matrix(0, n_ct, n_ct, dimnames = list(ct_names, ct_names))

for (i in seq_len(n_ct)) {
  for (j in seq_len(n_ct)) {
    overlap_mat[i, j] <- morisita_horn(ct_clone_mat[, i], ct_clone_mat[, j])
  }
}

overlap_df <- as.data.frame(overlap_mat) %>%
  tibble::rownames_to_column("cell_type_1") %>%
  tidyr::pivot_longer(-cell_type_1,
                      names_to = "cell_type_2", values_to = "overlap")

plot_overlap <- ggplot(overlap_df,
                       aes(x = cell_type_1, y = cell_type_2, fill = overlap)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(overlap, 2)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#d62728",
                       name = "Morisita-Horn") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL,
       title = "TCR clone overlap between T cell subtypes")

ggsave("clone_overlap_celltypes.pdf", plot_overlap, width = 7, height = 6)
write.csv(overlap_df, file.path(out_dir, "clone_overlap.csv"), row.names = FALSE)
