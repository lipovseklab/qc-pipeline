#!/usr/bin/env python3
"""Cluster comparison plots and clustree input for one topometry AnnData, for every clustering series it has:

    pca      pca_leiden_res*        PCA-based latent space
    topo     topo_clusters_res*     Topometry - spectral latent space
    topo_ms  topo_clusters_ms_res*  Topometry - multiscale latent space

Writes, in the current directory:
    cluster_comparison_on_embedding_<series>.png   the projection coloured by --color and by every resolution in the series
    cluster_table.csv                              cells x resolutions of all series, the input for bin/clustree_series.R
    clustree_series.txt                            one tab-separated line per series: name, column prefix, plot title
"""
import argparse
import sys
from pathlib import Path

import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parent))   # bin/ modules, wherever this runs from
from cluster_series import SERIES, plot_comparison, series_columns, write_clustree_input


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--adata', required=True, help='topometry AnnData (adata_topometry_*.h5ad)')
    p.add_argument('--color', default='ctype', help='obs column shown next to the clusterings (skipped if missing)')
    p.add_argument('--projection', default='projection', help='obsm key of the embedding to plot on')
    return p


def main():
    args = build_parser().parse_args()

    # only obs, obsm and uns are needed, so leave the expression matrix on disk
    adata = ad.read_h5ad(args.adata, backed='r')
    print(f'latent space chosen by topometry: {adata.uns["basis"]}')

    colours = []
    if args.color in adata.obs.columns:
        colours.append(args.color)
    else:
        print(f'WARNING: {args.color} is not in adata.obs; plotting the clusterings only')

    found = []
    for name, prefix, title in SERIES:
        cols = series_columns(adata.obs, prefix)
        if not cols:
            print(f'{name}: no {prefix}* columns, skipped')
            continue
        print(f'{name}: {len(cols)} resolutions ({prefix}*)')
        plot_comparison(adata, colours, cols, args.projection, f'cluster_comparison_on_embedding_{name}.png')
        found.append((name, prefix, title))

    if not found:
        raise SystemExit(f'None of the clustering series {[s[1] + "*" for s in SERIES]} are in {args.adata}')
    write_clustree_input(adata, found)


if __name__ == '__main__':
    main()
