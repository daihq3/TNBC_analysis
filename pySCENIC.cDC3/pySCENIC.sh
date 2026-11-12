python /home/data/t180330/miniconda3/envs/pyscenic/bin/arboreto_with_multiprocessing.py \
  --num_workers 20 \
  --output /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/grn.output.tsv \
  --method grnboost2 \
  /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/cell_type.loom \
  /home/data/t180330/cisTarget_HG38/allTFs_hg38.txt



pyscenic ctx \
  /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/grn.output.tsv \
  /home/data/t180330/cisTarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  /home/data/t180330/cisTarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --annotations_fname /home/data/t180330/cisTarget_HG38/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/cell_type.loom \
  --mode "custom_multiprocessing" \
  --output /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/reg.csv \
  --num_workers 60 \
  --mask_dropouts



pyscenic aucell \
  /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/cell_type.loom \
  /home/data/t180330/TNBC.ref.IL27.260115/Seurat_csv/reg.csv \
  --output SCENIC.csv \
  --num_workers 80

