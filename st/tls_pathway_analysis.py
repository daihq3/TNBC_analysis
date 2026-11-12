"""
TLS pathway analysis

Identify upregulated pathways in TLS-like regions using spatial transcriptomics.

Workflow:
  1. Define TLS spots from TLS_score
  2. Differential expression: TLS vs non-TLS spots (Wilcoxon rank-sum)
  3. Over-representation analysis (ORA) with hypergeometric test
  4. Visualization of top pathways

Author: Huiqi Dai
Date: 2026-08-08
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import gseapy as gp
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, hypergeom


def de_spots(adata, group_col='is_TLS', group_true=True):
    """
    Perform Wilcoxon rank-sum test for each gene between two spot groups.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object.
    group_col : str
        Column in adata.obs defining groups.
    group_true : bool
        Value in group_col corresponding to the test group (e.g., TLS).

    Returns
    -------
    pd.DataFrame
        Differential expression results with columns:
        gene, mean_group, mean_other, log2FC, pvalue, padj
    """
    group_mask = adata.obs[group_col] == group_true
    expr = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    results = []
    for i, gene in enumerate(adata.var_names):
        group_expr = expr[group_mask, i]
        other_expr = expr[~group_mask, i]

        # Skip low-expression or zero-variance genes
        if group_expr.std() == 0 or other_expr.std() == 0:
            continue
        if group_expr.mean() == 0 and other_expr.mean() == 0:
            continue

        stat, pval = mannwhitneyu(group_expr, other_expr, alternative='two-sided')

        results.append({
            'gene': gene,
            'mean_TLS': group_expr.mean(),
            'mean_nonTLS': other_expr.mean(),
            'log2FC': np.log2((group_expr.mean() + 1e-6) / (other_expr.mean() + 1e-6)),
            'pvalue': pval
        })

    results = pd.DataFrame(results)

    # FDR correction
    try:
        from statsmodels.stats.multitest import multipletests
        results['padj'] = multipletests(results['pvalue'], method='fdr_bh')[1]
    except ImportError:
        results['padj'] = results['pvalue']
        print("Warning: statsmodels not available; using raw p-values as padj.")

    results = results.sort_values('log2FC', ascending=False)
    return results


def read_gmt(gmt_file):
    """
    Read a GMT gene set file.

    Parameters
    ----------
    gmt_file : str
        Path to GMT file.

    Returns
    -------
    dict
        {pathway_name: set(genes)}
    """
    gmt_file = os.path.expanduser(gmt_file)
    pathways = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            pathway_name = parts[0]
            genes = parts[2:]
            pathways[pathway_name] = set(genes)
    return pathways


def ora_hypergeom(
    de_results,
    gmt_file,
    log2fc_thr=0.5,
    padj_thr=0.05,
    min_pathway_size=5
):
    """
    Over-representation analysis using hypergeometric test.

    Parameters
    ----------
    de_results : pd.DataFrame
        Output from de_spots(). Must contain 'gene', 'log2FC', 'padj'.
    gmt_file : str
        Path to GMT file.
    log2fc_thr : float
        log2 fold change threshold for up/down regulated genes.
    padj_thr : float
        Adjusted p-value threshold.
    min_pathway_size : int
        Minimum number of background genes in pathway to test.

    Returns
    -------
    pd.DataFrame
        ORA results with up/down enrichment p-values and FDR.
    """
    pathways = read_gmt(gmt_file)

    background_genes = set(de_results['gene'].values)

    up_genes = set(de_results[
        (de_results['log2FC'] > log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ]['gene'].values)

    down_genes = set(de_results[
        (de_results['log2FC'] < -log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ]['gene'].values)

    N = len(background_genes)
    K_up = len(up_genes)
    K_down = len(down_genes)

    results = []
    for pathway, pathway_genes in pathways.items():
        pathway_genes = pathway_genes & background_genes
        M = len(pathway_genes)

        if M < min_pathway_size:
            continue

        # Up-regulated enrichment
        x_up = len(pathway_genes & up_genes)
        p_up = hypergeom.sf(x_up - 1, N, M, K_up) if x_up > 0 else 1.0

        # Down-regulated enrichment
        x_down = len(pathway_genes & down_genes)
        p_down = hypergeom.sf(x_down - 1, N, M, K_down) if x_down > 0 else 1.0

        results.append({
            'pathway': pathway,
            'pathway_size': M,
            'up_genes_in_pathway': x_up,
            'down_genes_in_pathway': x_down,
            'up_pvalue': p_up,
            'down_pvalue': p_down,
            'up_genes': ','.join(sorted(pathway_genes & up_genes)),
            'down_genes': ','.join(sorted(pathway_genes & down_genes))
        })

    results = pd.DataFrame(results)

    # FDR correction
    try:
        from statsmodels.stats.multitest import multipletests
        results['up_padj'] = multipletests(results['up_pvalue'], method='fdr_bh')[1]
        results['down_padj'] = multipletests(results['down_pvalue'], method='fdr_bh')[1]
    except ImportError:
        results['up_padj'] = results['up_pvalue']
        results['down_padj'] = results['down_pvalue']
        print("Warning: statsmodels not available; using raw p-values as padj.")

    return results


def gseapy_ora(
    de_results,
    gmt_file,
    log2fc_thr=0.5,
    padj_thr=0.05,
    outdir=None,
    min_size=5,
    max_size=1000
):
    """
    Over-representation analysis using gseapy.enrich().

    Runs Fisher's exact test for up- and down-regulated gene sets
    against the GMT gene sets.

    Parameters
    ----------
    de_results : pd.DataFrame
        Output from de_spots().
    gmt_file : str
        Path to GMT file.
    log2fc_thr : float
        log2FC threshold.
    padj_thr : float
        Adjusted p-value threshold.
    outdir : str, optional
        Output directory for gseapy reports.
    min_size : int
        Minimum gene set size.
    max_size : int
        Maximum gene set size.

    Returns
    -------
    dict
        {'up': pd.DataFrame, 'down': pd.DataFrame, 'enr_up': Enrichr obj,
         'enr_down': Enrichr obj}
    """
    gmt_file = os.path.expanduser(gmt_file)
    background = de_results['gene'].tolist()

    up_genes = de_results[
        (de_results['log2FC'] > log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ]['gene'].tolist()

    down_genes = de_results[
        (de_results['log2FC'] < -log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ]['gene'].tolist()

    print(f"Up genes: {len(up_genes)}, Down genes: {len(down_genes)}")

    # Up-regulated ORA
    enr_up = gp.enrichr(
        gene_list=up_genes,
        gene_sets=gmt_file,
        background=background,
        outdir=f"{outdir}/up" if outdir else None,
        no_plot=True
    )
    df_up = enr_up.results.copy()
    df_up = df_up[df_up['Overlap'].apply(
        lambda x: int(str(x).split('/')[0]) >= min_size
    )].copy()

    # Down-regulated ORA
    enr_down = gp.enrichr(
        gene_list=down_genes,
        gene_sets=gmt_file,
        background=background,
        outdir=f"{outdir}/down" if outdir else None,
        no_plot=True
    )
    df_down = enr_down.results.copy()
    df_down = df_down[df_down['Overlap'].apply(
        lambda x: int(str(x).split('/')[0]) >= min_size
    )].copy()

    return {
        'up': df_up.sort_values('Adjusted P-value'),
        'down': df_down.sort_values('Adjusted P-value'),
        'enr_up': enr_up,
        'enr_down': enr_down
    }


def gseapy_prerank(
    de_results,
    gmt_file,
    outdir=None,
    min_size=5,
    max_size=1000,
    permutation_num=1000
):
    """
    Preranked GSEA using gseapy.prerank().

    Genes are ranked by log2FC (TLS vs non-TLS).

    Parameters
    ----------
    de_results : pd.DataFrame
        Output from de_spots().
    gmt_file : str
        Path to GMT file.
    outdir : str, optional
        Output directory for gseapy reports.
    min_size : int
        Minimum gene set size.
    max_size : int
        Maximum gene set size.
    permutation_num : int
        Number of permutations.

    Returns
    -------
    dict
        {'res': pd.DataFrame, 'prerank': Prerank obj}
    """
    gmt_file = os.path.expanduser(gmt_file)

    # Build ranked list: gene -> log2FC
    ranking = de_results[['gene', 'log2FC']].dropna()
    ranking = ranking.drop_duplicates(subset='gene')
    ranking = ranking.set_index('gene')['log2FC']

    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=gmt_file,
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        no_plot=True
    )

    df = pre_res.res2d.copy()
    df = df.sort_values('FDR q-val')

    return {
        'res': df,
        'prerank': pre_res
    }


def plot_gseapy_barplot(
    ora_df,
    direction='up',
    n_top=15,
    pval_col='Adjusted P-value',
    save_path=None
):
    """
    Barplot of top enriched pathways from gseapy ORA results.

    Parameters
    ----------
    ora_df : pd.DataFrame
        gseapy enrich results DataFrame.
    direction : str
        'up' or 'down'.
    n_top : int
        Number of top pathways to show.
    pval_col : str
        Column name for p-value.
    save_path : str, optional
        Path to save figure.
    """
    df = ora_df.head(n_top).copy()
    df['neg_log10_padj'] = -np.log10(df[pval_col].replace(0, 1e-300))

    color = '#C82F71' if direction == 'up' else '#7473A6'
    title = f'Top {n_top} {"up" if direction == "up" else "down"}regulated pathways in TLS'

    # Clean pathway names
    df['pathway_clean'] = df['Term'].str.replace('^HALLMARK_', '', regex=True).str.replace('_', ' ')

    plt.figure(figsize=(10, 0.5 * n_top + 2))
    plt.barh(df['pathway_clean'][::-1], df['neg_log10_padj'][::-1], color=color)
    plt.xlabel('-log10(FDR)')
    plt.ylabel('Pathway')
    plt.title(title, fontsize=13, fontweight='bold')
    plt.axvline(x=-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def get_top_pathways(ora_results, direction='up', n_top=20):
    """
    Get top enriched pathways.

    Parameters
    ----------
    ora_results : pd.DataFrame
        Output from ora_hypergeom().
    direction : str
        'up' or 'down'.
    n_top : int
        Number of top pathways to return.

    Returns
    -------
    pd.DataFrame
    """
    if direction == 'up':
        df = ora_results[
            (ora_results['up_genes_in_pathway'] > 0)
        ].sort_values('up_padj').head(n_top)
        return df[['pathway', 'pathway_size', 'up_genes_in_pathway',
                   'up_pvalue', 'up_padj', 'up_genes']]
    else:
        df = ora_results[
            (ora_results['down_genes_in_pathway'] > 0)
        ].sort_values('down_padj').head(n_top)
        return df[['pathway', 'pathway_size', 'down_genes_in_pathway',
                   'down_pvalue', 'down_padj', 'down_genes']]


def plot_top_pathways(ora_results, direction='up', n_top=15, save_path=None):
    """
    Barplot of top enriched pathways.
    """
    if direction == 'up':
        df = ora_results[
            (ora_results['up_genes_in_pathway'] > 0)
        ].sort_values('up_padj').head(n_top).copy()
        df['-log10_padj'] = -np.log10(df['up_padj'].replace(0, 1e-300))
        color = '#E74C3C'
        title = f'Top {n_top} upregulated pathways in TLS-like regions'
    else:
        df = ora_results[
            (ora_results['down_genes_in_pathway'] > 0)
        ].sort_values('down_padj').head(n_top).copy()
        df['-log10_padj'] = -np.log10(df['down_padj'].replace(0, 1e-300))
        color = '#3498DB'
        title = f'Top {n_top} downregulated pathways in TLS-like regions'

    # Clean pathway names
    df['pathway_clean'] = df['pathway'].str.replace('^HALLMARK_', '', regex=True).str.replace('_', ' ')

    plt.figure(figsize=(10, 0.5 * n_top + 2))
    plt.barh(df['pathway_clean'][::-1], df['-log10_padj'][::-1], color=color)
    plt.xlabel('-log10(FDR)')
    plt.ylabel('Pathway')
    plt.title(title)
    plt.axvline(x=-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def plot_gseapy_enrichment(
    prerank_obj,
    terms,
    figsize=(6, 4),
    save_prefix='GSEA'
):
    """
    Plot classic GSEA running-sum enrichment curves for selected pathways.

    Uses gseapy's built-in plot() method on the prerank object.

    Parameters
    ----------
    prerank_obj : gseapy.prerank.Prerank
        Output from gseapy_prerank()['prerank'].
    terms : list of str
        Pathway names (Term column) to plot.
    figsize : tuple
        Figure size per pathway.
    save_prefix : str
        Prefix for saved PDFs.
    """
    for term in terms:
        try:
            ax = prerank_obj.plot(
                term=term,
                figsize=figsize,
                ofname=f"{save_prefix}_{term.replace(' ', '_').replace('/', '_')}.pdf"
            )
            if ax is not None:
                plt.close()
                print(f"  Saved: {save_prefix}_{term}.pdf")
        except Exception as e:
            print(f"  Skipped {term}: {e}")


def plot_gseapy_dotplot(
    prerank_df,
    n_top=15,
    save_path=None
):
    """
    Dotplot of top enriched pathways from preranked GSEA.

    Shows NES (normalized enrichment score) colored by FDR q-value.
    Positive NES = enriched at top of ranked list (TLS-high).
    Negative NES = enriched at bottom (non-TLS-high).

    Parameters
    ----------
    prerank_df : pd.DataFrame
        Output from gseapy_prerank()['res'].
    n_top : int
        Number of top pathways (by absolute NES) to show.
    save_path : str, optional
        Path to save figure.
    """
    df = prerank_df.copy()

    # Select top by absolute NES
    df['abs_NES'] = df['NES'].abs()
    df_top = df.nlargest(n_top, 'abs_NES').copy()

    # Sort by NES for display
    df_top = df_top.sort_values('NES')
    df_top['-log10_FDR'] = -np.log10(df_top['FDR q-val'].replace(0, 1e-300))

    # Clean names
    df_top['pathway_clean'] = df_top['Term'].str.replace(
        '^HALLMARK_', '', regex=True
    ).str.replace('_', ' ')

    fig, ax = plt.subplots(figsize=(10, 0.4 * n_top + 2))

    scatter = ax.scatter(
        df_top['NES'], df_top['pathway_clean'],
        s=df_top['-log10_FDR'] * 15,
        c=df_top['NES'],
        cmap='RdBu_r',
        edgecolors='black', linewidths=0.5,
        vmin=-2.5, vmax=2.5
    )

    ax.axvline(x=0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_ylabel('')
    ax.set_title(f'Preranked GSEA: top {n_top} pathways', fontsize=13,
                 fontweight='bold')

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label('NES', fontsize=10)

    # Size legend
    sizes = [2, 5, 10]
    for s in sizes:
        ax.scatter([], [], s=s * 15, c='grey', alpha=0.5,
                   label=f'-log10(FDR) = {s}')
    ax.legend(title='Significance', loc='lower right', fontsize=8)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def plot_volcano(
    de_results,
    log2fc_thr=0.5,
    padj_thr=0.05,
    n_label=20,
    point_size=1.2,
    label_size=9,
    xlim=(-6, 6),
    save_path=None
):
    """
    Volcano plot for TLS vs non-TLS differential expression.

    Matches R ggplot2 style:
    - Up = #C82F71 (red), Down = #7473A6 (blue), NS = grey80
    - Filled label boxes colored by direction, white italic text
    - Black segment lines connecting labels to points
    - y-axis clamped at 170 to avoid Inf

    Parameters
    ----------
    de_results : pd.DataFrame
        Output from de_spots().
    log2fc_thr : float
        log2FC threshold.
    padj_thr : float
        Adjusted p-value threshold.
    n_label : int
        Number of top genes to label.
    point_size : float
        Scatter point size.
    label_size : float
        Font size for gene labels.
    xlim : tuple
        X-axis limits.
    save_path : str, optional
        Path to save figure.
    """
    col_up   = '#C82F71'
    col_down = '#7473A6'
    col_ns   = '#CCCCCC'  # grey80

    df = de_results.copy()
    df['neg_log10_p'] = -np.log10(df['padj'].replace(0, 1e-300))
    df['neg_log10_p'] = df['neg_log10_p'].clip(upper=170)

    # Classify direction
    df['direction'] = 'NS'
    df.loc[(df['log2FC'] > log2fc_thr) & (df['padj'] < padj_thr), 'direction'] = 'Up'
    df.loc[(df['log2FC'] < -log2fc_thr) & (df['padj'] < padj_thr), 'direction'] = 'Down'

    fig, ax = plt.subplots(figsize=(10, 8))

    # 1. Background NS points (grey, small)
    ns = df[df['direction'] == 'NS']
    ax.scatter(ns['log2FC'], ns['neg_log10_p'],
               c=col_ns, s=point_size * 8, alpha=0.6, edgecolors='none')

    # 2. Significant points
    up = df[df['direction'] == 'Up']
    dn = df[df['direction'] == 'Down']
    ax.scatter(up['log2FC'], up['neg_log10_p'],
               c=col_up, s=point_size * 10, alpha=0.7, edgecolors='none')
    ax.scatter(dn['log2FC'], dn['neg_log10_p'],
               c=col_down, s=point_size * 10, alpha=0.7, edgecolors='none')

    # 3. Top genes to label (by combined score)
    sig = df[df['direction'] != 'NS'].copy()
    sig['score'] = sig['log2FC'].abs() * sig['neg_log10_p']
    lab_df = sig.nlargest(n_label, 'score')

    # Larger labeled points on top
    color_map = {'Up': col_up, 'Down': col_down}
    for _, row in lab_df.iterrows():
        ax.scatter(row['log2FC'], row['neg_log10_p'],
                   c=color_map.get(row['direction'], col_ns),
                   s=40, alpha=1.0, edgecolors='black', linewidths=0.3, zorder=5)

    # 4. Filled label boxes with white italic text + black segments
    from adjustText import adjust_text
    texts = []
    for _, row in lab_df.iterrows():
        box_color = col_up if row['direction'] == 'Up' else col_down
        texts.append(ax.annotate(
            row['gene'],
            (row['log2FC'], row['neg_log10_p']),
            textcoords='offset points',
            xytext=(12, 12),
            fontsize=label_size,
            fontstyle='italic',
            color='white',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', fc=box_color, ec='none', alpha=0.9),
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
            zorder=6
        ))

    # Adjust text positions to reduce overlap
    try:
        adjust_text(texts, ax=ax, force=(2, 2),
                    arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    except Exception:
        pass  # labels still placed if adjustText fails

    # 5. Threshold lines
    ax.axhline(y=-np.log10(padj_thr), color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(x=log2fc_thr, color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(x=-log2fc_thr, color='grey', linestyle='--', linewidth=0.8)

    # 6. Axes and theme
    ax.set_xlim(xlim)
    ax.set_xlabel(r'$\log_2$(fold change)', fontsize=13)
    ax.set_ylabel(r'$-\log_{10}$(P value)', fontsize=13)
    ax.set_title('TLS vs non-TLS spots', fontsize=13, fontweight='bold', loc='center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=11)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def plot_deg_heatmap(
    adata,
    de_results,
    log2fc_thr=0.5,
    padj_thr=0.05,
    n_top=25,
    tls_col='is_TLS',
    sample_col='sample_id',
    max_spots=2000,
    save_path=None
):
    """
    Heatmap of top DEGs: TLS vs non-TLS spots.

    Shows top up- and down-regulated genes (z-scored expression)
    with spots sorted by TLS status and colored by sample.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object.
    de_results : pd.DataFrame
        Output from de_spots().
    log2fc_thr : float
        log2FC threshold.
    padj_thr : float
        Adjusted p-value threshold.
    n_top : int
        Number of top up + top down genes to show.
    tls_col : str
        Column for TLS label.
    sample_col : str
        Column for sample ID (used for column color bar).
    max_spots : int
        Maximum number of spots to subsample for plotting.
    save_path : str, optional
        Path to save figure.
    """
    from scipy.stats import zscore
    import seaborn as sns

    # 1. Select top genes
    sig_up = de_results[
        (de_results['log2FC'] > log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ].nlargest(n_top, 'log2FC')

    sig_down = de_results[
        (de_results['log2FC'] < -log2fc_thr) &
        (de_results['padj'] < padj_thr)
    ].nsmallest(n_top, 'log2FC')

    top_genes = list(sig_up['gene']) + list(sig_down['gene'])
    top_genes = [g for g in top_genes if g in adata.var_names]
    # Deduplicate while preserving order
    seen = set()
    top_genes = [g for g in top_genes if not (g in seen or seen.add(g))]

    print(f"Heatmap genes: {len(top_genes)} ({len(sig_up)} up + {len(sig_down)} down)")

    # 2. Extract expression and z-score per gene
    expr = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    gene_matrix = np.column_stack([expr[:, var_idx[g]] for g in top_genes])
    gene_matrix_z = zscore(gene_matrix, axis=0, nan_policy='omit')
    gene_matrix_z = np.nan_to_num(gene_matrix_z, nan=0.0)

    # 3. Sort spots: non-TLS first, then TLS
    tls_labels = adata.obs[tls_col].values.astype(bool)
    sort_idx = np.argsort(tls_labels)  # False (non-TLS) first

    # Subsample if too many spots
    if len(sort_idx) > max_spots:
        rng = np.random.RandomState(42)
        tls_idx = sort_idx[tls_labels[sort_idx]]
        nontls_idx = sort_idx[~tls_labels[sort_idx]]
        n_each = max_spots // 2
        if len(tls_idx) > n_each:
            tls_idx = rng.choice(tls_idx, n_each, replace=False)
        if len(nontls_idx) > n_each:
            nontls_idx = rng.choice(nontls_idx, n_each, replace=False)
        sort_idx = np.concatenate([nontls_idx, tls_idx])
        sort_idx = np.sort(sort_idx)  # keep original order within group
        # Re-sort by TLS within subsample
        sub_tls = tls_labels[sort_idx]
        sort_idx = sort_idx[np.argsort(sub_tls)]

    # 4. Build expression matrix for plotting
    mat_z = gene_matrix_z[sort_idx].T  # genes x spots
    tls_sorted = tls_labels[sort_idx]

    # 5. Plot with matplotlib
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch
    import matplotlib.gridspec as gridspec

    n_genes = len(top_genes)
    n_spots = len(sort_idx)
    fig_height = max(6, n_genes * 0.35)
    fig_width = max(10, n_spots * 0.005 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(
        2, 1, height_ratios=[0.3, n_genes * 0.35],
        hspace=0.02
    )

    # Top bar: TLS status color strip
    ax_bar = fig.add_subplot(gs[0])
    bar_colors = np.array([
        [0.784, 0.184, 0.443] if t else [0.455, 0.451, 0.651]
        for t in tls_sorted
    ])  # #C82F71 vs #7473A6
    ax_bar.imshow(
        bar_colors[np.newaxis, :, :],
        aspect='auto', interpolation='none'
    )
    ax_bar.set_yticks([0])
    ax_bar.set_yticklabels(['TLS'])
    ax_bar.set_xticks([])
    ax_bar.tick_params(axis='y', labelsize=9)

    # Main heatmap
    ax_heat = fig.add_subplot(gs[1])
    cmap = plt.cm.RdBu_r
    im = ax_heat.imshow(
        mat_z, aspect='auto', cmap=cmap,
        vmin=-2.5, vmax=2.5, interpolation='none'
    )
    ax_heat.set_yticks(range(n_genes))
    ax_heat.set_yticklabels(top_genes, fontsize=9)
    ax_heat.set_xticks([])
    ax_heat.set_xlabel(f'Spots (n={n_spots})', fontsize=11)
    ax_heat.set_title('Top DEGs: TLS vs non-TLS', fontsize=13,
                      fontweight='bold', pad=10)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.5, pad=0.02)
    cbar.set_label('Z-score', fontsize=10)

    # Legend
    legend_elements = [
        Patch(facecolor='#C82F71', label='TLS'),
        Patch(facecolor='#7473A6', label='non-TLS')
    ]
    ax_bar.legend(
        handles=legend_elements, loc='upper right',
        bbox_to_anchor=(1.12, 1.0), fontsize=9, frameon=True
    )

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    plt.show()
    return top_genes


def run_tls_pathway_analysis(
    slide,
    gmt_file,
    tls_score_col='TLS_score',
    tls_quantile=0.90,
    log2fc_thr=0.5,
    padj_thr=0.05,
    save_prefix='TLS'
):
    """
    One-step wrapper: define TLS spots, DE, ORA, and save results.

    Parameters
    ----------
    slide : AnnData
        Spatial transcriptomics object with TLS_score in .obs.
    gmt_file : str
        Path to GMT file.
    tls_score_col : str
        Column name for TLS score.
    tls_quantile : float
        Quantile threshold to define TLS spots.
    log2fc_thr : float
        log2FC threshold for ORA.
    padj_thr : float
        Adjusted p-value threshold for ORA.
    save_prefix : str
        Prefix for output CSV files.

    Returns
    -------
    dict
        {'de_results': ..., 'ora_results': ..., 'top_up': ..., 'top_down': ...}
    """
    # 1. Define TLS spots
    threshold = slide.obs[tls_score_col].quantile(tls_quantile)
    slide.obs['is_TLS'] = slide.obs[tls_score_col] > threshold
    print(f"TLS spots: {slide.obs['is_TLS'].sum()} / {slide.n_obs}")

    # 2. Differential expression
    print("Running differential expression...")
    de_results = de_spots(slide, group_col='is_TLS', group_true=True)
    de_results.to_csv(f"{save_prefix}_DE_results.csv", index=False)

    # 3. gseapy ORA
    print("Running gseapy ORA...")
    gsea_ora = gseapy_ora(
        de_results, gmt_file,
        log2fc_thr=log2fc_thr, padj_thr=padj_thr,
        outdir=f"{save_prefix}_gseapy_ORA"
    )
    gsea_ora['up'].to_csv(f"{save_prefix}_ORA_up.csv", index=False)
    gsea_ora['down'].to_csv(f"{save_prefix}_ORA_down.csv", index=False)

    # 4. gseapy preranked GSEA
    print("Running gseapy preranked GSEA...")
    gsea_pr = gseapy_prerank(
        de_results, gmt_file,
        outdir=f"{save_prefix}_gseapy_GSEA"
    )
    gsea_pr['res'].to_csv(f"{save_prefix}_GSEA_prerank.csv", index=False)

    # 5. Plot ORA barplots
    plot_gseapy_barplot(gsea_ora['up'], direction='up', n_top=15,
                        save_path=f"{save_prefix}_ORA_up_pathways.pdf")
    plot_gseapy_barplot(gsea_ora['down'], direction='down', n_top=15,
                        save_path=f"{save_prefix}_ORA_down_pathways.pdf")

    # 5b. GSEA preranked dotplot
    plot_gseapy_dotplot(gsea_pr['res'], n_top=15,
                        save_path=f"{save_prefix}_GSEA_dotplot.pdf")

    # 5c. GSEA enrichment curves for top 5 pathways
    top_terms = gsea_pr['res'].nlargest(5, 'NES.abs')['Term'].tolist() if 'NES.abs' in gsea_pr['res'].columns else \
                gsea_pr['res'].assign(a=gsea_pr['res']['NES'].abs()).nlargest(5, 'a')['Term'].tolist()
    print(f"Plotting enrichment curves for top 5 pathways...")
    plot_gseapy_enrichment(gsea_pr['prerank'], top_terms,
                           save_prefix=f"{save_prefix}_GSEA_curve")

    # 6. Volcano plot
    plot_volcano(de_results, log2fc_thr=log2fc_thr, padj_thr=padj_thr,
                 save_path=f"{save_prefix}_volcano.pdf")

    # 7. DEG heatmap
    plot_deg_heatmap(slide, de_results, log2fc_thr=log2fc_thr,
                     padj_thr=padj_thr, n_top=25,
                     save_path=f"{save_prefix}_DEG_heatmap.pdf")

    return {
        'de_results': de_results,
        'gseapy_ora': gsea_ora,
        'gseapy_prerank': gsea_pr
    }


# ==============================================================================
# Multi-slide merge and merged TLS analysis
# ==============================================================================


def merge_tls_slides(slides_dict, celltype_cols=None):
    """
    Merge multiple Visium slides into a single AnnData.

    Only genes present in all slides (intersection) are kept.
    Optionally copies cell2location-derived cell-type abundance columns
    from .obs into the merged object.

    Parameters
    ----------
    slides_dict : dict
        {slide_name: AnnData}
    celltype_cols : list of str, optional
        Columns to copy from .obs. If None, all .obs columns are kept.

    Returns
    -------
    anndata.AnnData
        Merged AnnData with 'sample_id' column in .obs.
    """
    merged_list = []
    for name, adata in slides_dict.items():
        adata_copy = adata.copy()
        adata_copy.obs['sample_id'] = name

        if celltype_cols is not None:
            keep_cols = [c for c in celltype_cols if c in adata_copy.obs.columns]
            adata_copy.obs = adata_copy.obs[['sample_id'] + keep_cols]

        merged_list.append(adata_copy)

    merged = ad.concat(merged_list, join='inner', merge='same', label='sample_id')
    merged.obs_names_make_unique()
    print(f"Merged {len(slides_dict)} slides: {merged.n_obs} spots, {merged.n_vars} genes")
    return merged


def define_tls_per_slide(adata, tls_score_col='TLS_score', tls_quantile=0.90,
                         output_col='is_TLS'):
    """
    Define TLS spots per-slide using per-slide quantile threshold.

    Each slide gets its own TLS threshold so that the top N% of spots
    in every slide are labeled TLS.

    Parameters
    ----------
    adata : AnnData
        Merged AnnData with 'sample_id' and TLS_score column.
    tls_score_col : str
        Column name for TLS score.
    tls_quantile : float
        Quantile threshold (applied per slide).
    output_col : str
        Name of the output boolean column.

    Returns
    -------
    anndata.AnnData
        AnnData with `output_col` in .obs.
    """
    is_tls = np.zeros(adata.n_obs, dtype=bool)

    for slide_name in adata.obs['sample_id'].unique():
        mask = adata.obs['sample_id'] == slide_name
        slide_scores = adata.obs.loc[mask, tls_score_col]
        threshold = slide_scores.quantile(tls_quantile)
        is_tls[mask] = slide_scores > threshold
        n_tls = (is_tls[mask]).sum()
        print(f"  {slide_name}: {n_tls} TLS spots / {mask.sum()} total "
              f"(threshold={threshold:.3f})")

    adata.obs[output_col] = is_tls
    return adata


def run_merged_tls_pathway_analysis(
    slides_dict,
    gmt_file,
    tls_score_col='TLS_score',
    tls_quantile=0.90,
    celltype_cols=None,
    log2fc_thr=0.5,
    padj_thr=0.05,
    save_prefix='merged_TLS'
):
    """
    Merge slides, define TLS per-slide, then run DE + ORA on combined data.

    Parameters
    ----------
    slides_dict : dict
        {slide_name: AnnData}
    gmt_file : str
        Path to GMT gene set file.
    tls_score_col : str
        Column name for TLS score.
    tls_quantile : float
        Per-slide quantile threshold for TLS definition.
    celltype_cols : list of str, optional
        Cell-type columns to carry through merge.
    log2fc_thr : float
        log2FC threshold for ORA.
    padj_thr : float
        Adjusted p-value threshold for ORA.
    save_prefix : str
        Prefix for output CSV and PDF files.

    Returns
    -------
    dict
        {'merged_adata': ..., 'de_results': ...,
         'ora_results': ..., 'top_up': ..., 'top_down': ...}
    """
    # 1. Merge slides
    print("Merging slides...")
    merged = merge_tls_slides(slides_dict, celltype_cols=celltype_cols)

    # 2. Define TLS per slide
    print("Defining TLS spots per slide...")
    merged = define_tls_per_slide(merged, tls_score_col=tls_score_col,
                                   tls_quantile=tls_quantile)
    total_tls = merged.obs['is_TLS'].sum()
    print(f"Total TLS spots: {total_tls} / {merged.n_obs}")

    # 3. Differential expression
    print("Running differential expression (TLS vs non-TLS)...")
    de_results = de_spots(merged, group_col='is_TLS', group_true=True)
    de_results.to_csv(f"{save_prefix}_DE_results.csv", index=False)

    # 4. gseapy ORA
    print("Running gseapy ORA...")
    gsea_ora = gseapy_ora(
        de_results, gmt_file,
        log2fc_thr=log2fc_thr, padj_thr=padj_thr,
        outdir=f"{save_prefix}_gseapy_ORA"
    )
    gsea_ora['up'].to_csv(f"{save_prefix}_ORA_up.csv", index=False)
    gsea_ora['down'].to_csv(f"{save_prefix}_ORA_down.csv", index=False)

    # 5. gseapy preranked GSEA
    print("Running gseapy preranked GSEA...")
    gsea_pr = gseapy_prerank(
        de_results, gmt_file,
        outdir=f"{save_prefix}_gseapy_GSEA"
    )
    gsea_pr['res'].to_csv(f"{save_prefix}_GSEA_prerank.csv", index=False)

    # 6. Plot ORA barplots
    plot_gseapy_barplot(gsea_ora['up'], direction='up', n_top=15,
                        save_path=f"{save_prefix}_ORA_up_pathways.pdf")
    plot_gseapy_barplot(gsea_ora['down'], direction='down', n_top=15,
                        save_path=f"{save_prefix}_ORA_down_pathways.pdf")

    # 6b. GSEA preranked dotplot
    plot_gseapy_dotplot(gsea_pr['res'], n_top=15,
                        save_path=f"{save_prefix}_GSEA_dotplot.pdf")

    # 6c. GSEA enrichment curves for top 5 pathways
    top_terms = gsea_pr['res'].nlargest(5, 'NES.abs')['Term'].tolist() if 'NES.abs' in gsea_pr['res'].columns else \
                gsea_pr['res'].assign(a=gsea_pr['res']['NES'].abs()).nlargest(5, 'a')['Term'].tolist()
    print(f"Plotting enrichment curves for top 5 pathways...")
    plot_gseapy_enrichment(gsea_pr['prerank'], top_terms,
                           save_prefix=f"{save_prefix}_GSEA_curve")

    # 7. Volcano plot
    plot_volcano(de_results, log2fc_thr=log2fc_thr, padj_thr=padj_thr,
                 save_path=f"{save_prefix}_volcano.pdf")

    # 8. DEG heatmap
    plot_deg_heatmap(merged, de_results, log2fc_thr=log2fc_thr,
                     padj_thr=padj_thr, n_top=25,
                     save_path=f"{save_prefix}_DEG_heatmap.pdf")

    return {
        'merged_adata': merged,
        'de_results': de_results,
        'gseapy_ora': gsea_ora,
        'gseapy_prerank': gsea_pr
    }


# ==============================================================================
# Example usage
# ==============================================================================

# import os
# os.chdir('data_transfer/TLS')
#
# from tls_pathway_analysis import (
#     run_tls_pathway_analysis,
#     run_merged_tls_pathway_analysis
# )
#
# gmt_file = "~/202311_tnbc_caf_05_20240521/202311_tnbc_caf_05_20240812/cell_type_DEG/T/pathway/h.all.v2024.1.Hs.symbols.gmt"
#
# # --- Single-slide analysis ---
# results = run_tls_pathway_analysis(
#     slide_4465,
#     gmt_file=gmt_file,
#     tls_score_col='TLS_score',
#     save_prefix='CID4465_TLS'
# )
#
# # --- Multi-slide merged analysis ---
# slides = {
#     'CID4465':           slide_4465,
#     'V1_blockA_section2': slide_v1,
#     'Visium_FFPE':       slide_ffpe
# }
#
# results_merged = run_merged_tls_pathway_analysis(
#     slides,
#     gmt_file=gmt_file,
#     tls_score_col='TLS_score',
#     tls_quantile=0.90,
#     save_prefix='merged_TLS'
# )
#
# print(results_merged['gseapy_ora']['up'].head(10))
# print(results_merged['gseapy_prerank']['res'].head(10))
