"""
Radial distance analysis from TLS-like regions.

Python implementation inspired by semla::RadialDistance:
- Distance = 0 represents spots immediately exterior to TLS boundary.
- Positive distances increase outward from TLS.
- Distances > max_dist (default 1000 um) are discarded.
- Cell-type density is binned by radial distance and smoothed.

Author: Qoder
Date: 2026-08-08
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import UnivariateSpline


def sqrt_scale_distance(distances):
    """
    Apply signed square-root transformation to radial distances.

    d_sqrt = sign(d) * sqrt(|d|)

    This compresses large distances while preserving the inside/outside
    direction relative to the TLS boundary.
    """
    return np.sign(distances) * np.sqrt(np.abs(distances))


def compute_radial_distance_from_tls(
    slide,
    tls_col='is_TLS',
    max_dist=1000,
    k_neighbors=6,
    inside_negative=True,
    sqrt_scale=False,
    output_col='radial_distance_tls'
):
    """
    Compute radial distance from TLS boundary for each spot.

    By default, uses Convention 2: negative inside TLS, positive outside TLS.
    Distance 0 represents the TLS boundary.

    Parameters
    ----------
    slide : AnnData
        Spatial transcriptomics object.
    tls_col : str
        Column with TLS boolean label.
    max_dist : float
        Maximum distance to keep; spots beyond this are set to NaN.
    k_neighbors : int
        Number of nearest neighbors to determine TLS boundary.
    inside_negative : bool
        If True, interior TLS spots have negative distances (Convention 2).
        If False, interior TLS spots are set to 0 (exterior-only convention).
    sqrt_scale : bool
        If True, apply signed square-root transformation to distances.
    output_col : str
        Name of the output column in slide.obs.

    Returns
    -------
    AnnData
        Input object with added radial_distance_tls column.
    """
    coords = slide.obsm['spatial']
    tls_mask = slide.obs[tls_col].values.astype(bool)

    # Find TLS boundary spots: TLS spots with at least one non-TLS neighbor
    tree = cKDTree(coords)
    _, neighbors = tree.query(coords, k=k_neighbors)

    boundary_mask = np.zeros(len(coords), dtype=bool)
    for i in range(len(coords)):
        if tls_mask[i]:
            if np.any(~tls_mask[neighbors[i]]):
                boundary_mask[i] = True

    # Distance from each spot to nearest TLS boundary spot
    boundary_coords = coords[boundary_mask]
    if len(boundary_coords) == 0:
        slide.obs[output_col] = np.nan
        return slide

    boundary_tree = cKDTree(boundary_coords)
    distances, _ = boundary_tree.query(coords)

    # Convention 2: negative inside TLS, positive outside TLS
    radial_dist = distances.copy()
    if inside_negative:
        radial_dist[tls_mask] = -radial_dist[tls_mask]
    else:
        # Convention 1: exterior only, interior spots set to 0
        radial_dist[tls_mask] = 0.0

    # Discard distant spots (apply to absolute distance)
    radial_dist[np.abs(radial_dist) > max_dist] = np.nan

    # Optional signed square-root transformation
    if sqrt_scale:
        radial_dist = sqrt_scale_distance(radial_dist)

    slide.obs[output_col] = radial_dist
    return slide


def cell_density_by_distance(
    slide,
    celltype_col,
    distance_col='radial_distance_tls',
    bin_width=50,
    max_dist=1000
):
    """
    Compute mean cell-type abundance in radial distance bins.

    Parameters
    ----------
    slide : AnnData
        Spatial transcriptomics object.
    celltype_col : str
        Column in slide.obs with cell-type abundance.
    distance_col : str
        Column in slide.obs with radial distance.
    bin_width : float
        Width of radial distance bins.
    max_dist : float
        Maximum radial distance.

    Returns
    -------
    pd.DataFrame
        Binned density with columns: distance, abundance, celltype, n_spots.
    """
    df = pd.DataFrame({
        'distance': slide.obs[distance_col].values,
        'abundance': slide.obs[celltype_col].values
    })

    df = df.dropna()
    if len(df) == 0:
        return pd.DataFrame(columns=['distance', 'abundance', 'celltype', 'n_spots'])

    # Support negative distances (inside TLS) and positive distances (outside TLS)
    min_dist = np.floor(df['distance'].min() / bin_width) * bin_width
    max_dist_ceil = np.ceil(df['distance'].max() / bin_width) * bin_width
    bins = np.arange(min_dist, max_dist_ceil + bin_width, bin_width)
    df['distance_bin'] = pd.cut(df['distance'], bins=bins, include_lowest=True)

    binned = df.groupby('distance_bin').agg(
        abundance=('abundance', 'mean'),
        distance=('distance', 'mean'),
        n_spots=('abundance', 'size')
    ).reset_index()

    binned = binned.dropna()
    binned['celltype'] = celltype_col

    return binned


def smooth_radial_density(x, y, span=0.45, n_points=200):
    """
    Apply smoothing spline to radial density data (LOESS-like).
    """
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]

    if len(x) < 4:
        return x, y

    # Smoothing factor scales with number of points
    s = span * len(x)
    try:
        spline = UnivariateSpline(x, y, s=s)
        x_smooth = np.linspace(x.min(), x.max(), n_points)
        y_smooth = spline(x_smooth)
        return x_smooth, y_smooth
    except Exception as e:
        print(f"Smoothing failed: {e}; returning raw data")
        return x, y


def plot_radial_density(
    density_df,
    celltype,
    slide_name="sample",
    span=0.45,
    show_points=True,
    save_path=None
):
    """
    Plot cell-type density along radial distance with smoothing.
    """
    x = density_df['distance'].values
    y = density_df['abundance'].values

    x_smooth, y_smooth = smooth_radial_density(x, y, span=span)

    fig, ax = plt.subplots(figsize=(8, 5))

    if show_points:
        ax.scatter(x, y, alpha=0.5, s=30, color='gray', label='Binned mean')

    ax.plot(x_smooth, y_smooth, color='#E74C3C', linewidth=2, label='LOESS-like fit')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, label='TLS boundary')

    ax.set_xlabel('Radial distance from TLS (μm)')
    ax.set_ylabel(f'{celltype} abundance')
    ax.set_title(f'{celltype} density around TLS: {slide_name}')
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_multi_radial_density(
    slide,
    celltype_cols,
    slide_name="sample",
    bin_width=50,
    max_dist=1000,
    span=0.45,
    save_path=None
):
    """
    Plot multiple cell-type radial densities on the same axes.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.tab10(np.linspace(0, 1, len(celltype_cols)))

    for idx, ct in enumerate(celltype_cols):
        if ct not in slide.obs.columns:
            print(f"  Skipping {ct}: not found")
            continue

        density = cell_density_by_distance(
            slide, ct, bin_width=bin_width, max_dist=max_dist
        )
        if len(density) == 0:
            continue

        x = density['distance'].values
        y = density['abundance'].values
        x_smooth, y_smooth = smooth_radial_density(x, y, span=span)

        ax.plot(x_smooth, y_smooth, linewidth=2, color=colors[idx], label=ct)

    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
    ax.set_xlabel('Radial distance from TLS (μm)')
    ax.set_ylabel('Cell abundance')
    ax.set_title(f'Cell-type density around TLS: {slide_name}')
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def run_radial_distance_analysis(
    slide,
    celltype_cols,
    tls_col='is_TLS',
    slide_name="sample",
    bin_width=50,
    max_dist=1000,
    inside_negative=True,
    sqrt_scale=False,
    span=0.45,
    save_prefix='radial'
):
    """
    One-step wrapper: compute radial distance and plot densities.
    """
    # Ensure TLS label exists
    if tls_col not in slide.obs.columns:
        raise ValueError(
            f"'{tls_col}' not found in slide.obs. "
            "Please define TLS regions first (e.g., from TLS_score)."
        )

    # Compute radial distance
    slide = compute_radial_distance_from_tls(
        slide,
        tls_col=tls_col,
        max_dist=max_dist,
        inside_negative=inside_negative,
        sqrt_scale=sqrt_scale
    )

    all_density = []
    for ct in celltype_cols:
        if ct not in slide.obs.columns:
            print(f"  Skipping {ct}: not found in slide.obs")
            continue

        density = cell_density_by_distance(
            slide, ct, bin_width=bin_width, max_dist=max_dist
        )
        density.to_csv(f"{save_prefix}_{ct}_radial_density.csv", index=False)
        all_density.append(density)

        # Individual plot
        plot_radial_density(
            density,
            ct,
            slide_name=slide_name,
            span=span,
            save_path=f"{save_prefix}_{ct}_radial_density.pdf"
        )

    # Combined plot
    plot_multi_radial_density(
        slide,
        celltype_cols,
        slide_name=slide_name,
        bin_width=bin_width,
        max_dist=max_dist,
        span=span,
        save_path=f"{save_prefix}_multi_celltype_radial_density.pdf"
    )

    combined_density = pd.concat(all_density, ignore_index=True)
    combined_density.to_csv(f"{save_prefix}_combined_radial_density.csv", index=False)

    return slide, combined_density


# ==============================================================================
# Example usage
# ==============================================================================

# from tls_radial_distance import run_radial_distance_analysis
#
# # Make sure is_TLS exists first
# slide_4465.obs['is_TLS'] = slide_4465.obs['TLS_score'] > slide_4465.obs['TLS_score'].quantile(0.90)
#
# # Standard radial distance
# slide_4465, combined_density = run_radial_distance_analysis(
#     slide_4465,
#     celltype_cols=['B.cycling', 'Bmem', 'Bn.TCL1A', 'CD8.CXCL13', 'cDC3.LAMP3'],
#     slide_name='CID4465',
#     save_prefix='CID4465'
# )
#
# # Signed square-root scaled radial distance
# slide_4465, combined_density_sqrt = run_radial_distance_analysis(
#     slide_4465,
#     celltype_cols=['B.cycling', 'Bmem', 'Bn.TCL1A', 'CD8.CXCL13', 'cDC3.LAMP3'],
#     slide_name='CID4465',
#     sqrt_scale=True,
#     save_prefix='CID4465_sqrt'
# )
