"""Clustering series in a TopOMetry AnnData: which columns belong to which latent space, the cluster comparison plot,
and the input for the clustree plot (bin/clustree_series.R). Used by bin/combine_clusters.py (QC, outliers) and
bin/cluster_compare.py (cNMF).

Each series is one latent space clustered at several resolutions (bin/topometry.py):

    pca      pca_leiden_res<r>        PCA-based latent space
    topo     topo_clusters_res<r>     Topometry - spectral latent space
    topo_ms  topo_clusters_ms_res<r>  Topometry - multiscale latent space
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

# name, column prefix (followed by the resolution), title of the clustree plot
SERIES = [
    ('pca', 'pca_leiden_res', 'PCA-based latent space'),
    ('topo', 'topo_clusters_res', 'Topometry - spectral latent space'),
    ('topo_ms', 'topo_clusters_ms_res', 'Topometry - multiscale latent space'),
]


def series_columns(obs, prefix):
    """The columns of one series in order of resolution (e.g. topo_clusters_res0.2, ..., topo_clusters_res1.2)."""
    cols = [c for c in obs.columns if c.startswith(prefix)]
    resolution = {}
    for c in cols:
        try:
            resolution[c] = float(c[len(prefix):])
        except ValueError:                      # e.g. topo_clusters_res... when looking for pca_leiden_res: not ours
            pass
    return sorted(resolution, key=resolution.get)


def series_for_basis(basis):
    """The series TopOMetry's chosen latent space (adata.uns['basis']) belongs to.

    Multiscale bases contain 'ms' (e.g. 'ms_spectral_scaffold'); other spectral/TopOMetry bases ('spectral_scaffold',
    'topo...') are the spectral series; anything else is PCA.
    """
    basis = str(basis)
    if 'ms' in basis:
        name = 'topo_ms'
    elif 'spectral' in basis or 'topo' in basis:
        name = 'topo'
    else:
        name = 'pca'
    return next(s for s in SERIES if s[0] == name)


def plot_comparison(adata, colours, cols, projection, out):
    """The projection coloured by `colours` (e.g. cell types) and by every clustering in `cols`, saved to `out`.

    Clusterings that aren't categorical yet are shown as categories (not as a colour scale); adata.obs is left as it
    was, so a caller can still write the AnnData unchanged.
    """
    original = {c: adata.obs[c] for c in cols if not isinstance(adata.obs[c].dtype, pd.CategoricalDtype)}
    for c in original:
        adata.obs[c] = adata.obs[c].astype(str).astype('category')
    try:
        sc.pl.embedding(adata, basis=projection, color=colours + cols, legend_loc='on data',
                        legend_fontsize='xx-small', frameon=False, show=False)
        plt.savefig(out)
        plt.close()
    finally:
        for c, values in original.items():
            adata.obs[c] = values


def write_clustree_input(adata, series, table='cluster_table.csv', listing='clustree_series.txt'):
    """Write the cluster table (cells x resolutions of the given series) and the series list for clustree_series.R.

    `series` is a list of (name, prefix, title) tuples, each with at least one column in adata.obs.
    """
    cols = [c for _, prefix, _ in series for c in series_columns(adata.obs, prefix)]
    adata.obs[cols].astype(str).to_csv(table)
    with open(listing, 'w') as f:
        f.writelines(f'{name}\t{prefix}\t{title}\n' for name, prefix, title in series)
