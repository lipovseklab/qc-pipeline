#!/usr/bin/env python3
"""Add the cell types (and scSHC clusters) to the TopOMetry AnnData, plot the clusterings of TopOMetry's chosen latent
space next to them, and write the input for the clustree plot (bin/clustree_series.R).

Writes adata_clusters.h5ad, cluster_comparison_on_embedding.png, cluster_table.csv and clustree_series.txt.
"""
import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))   # bin/ modules, wherever this runs from
from cluster_series import plot_comparison, series_columns, series_for_basis, write_clustree_input


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--celltype-csv', required=True)
    p.add_argument('--scshc-csv', required=True)
    p.add_argument('--descriptor', required=False)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)
    ctype = pd.read_csv(args.celltype_csv, index_col=0, header=0)

    # rename old clusters & cell typing if they exist
    tag = 'qc' if args.descriptor == 'postqc' else f'before_{args.descriptor}'
    for key in ['shc_clusters', 'predicted_celltype']:
        if key in adata.obs:
            adata.obs[f'{key}_{tag}'] = adata.obs[key]
            del adata.obs[key]

    if os.path.exists(args.scshc_csv):
        shc = pd.read_csv(args.scshc_csv, index_col=0, header=0)
        adata.obs = adata.obs.join([ctype.predicted_celltype, shc.shc_clusters])
        adata.obs['shc_clusters'] = adata.obs['shc_clusters'].astype(int).astype('category')
    else:
        adata.obs = adata.obs.join(ctype.predicted_celltype)
    adata.obs['predicted_celltype'] = adata.obs['predicted_celltype'].astype('category')

    # the clusterings of the latent space TopOMetry chose (adata.uns['basis'])
    series = series_for_basis(adata.uns['basis'])
    cluster_cols = series_columns(adata.obs, series[1])
    print(f"latent space chosen by topometry: {adata.uns['basis']} -> {series[0]} ({len(cluster_cols)} resolutions)")
    if not cluster_cols:
        raise SystemExit(f"No {series[1]}* columns in {args.adata} for latent space {adata.uns['basis']}")

    colours = ['predicted_celltype'] + (['shc_clusters'] if 'shc_clusters' in adata.obs.columns else [])
    plot_comparison(adata, colours, cluster_cols, 'projection', 'cluster_comparison_on_embedding.png')
    write_clustree_input(adata, [series])

    adata.write('adata_clusters.h5ad')


if __name__ == '__main__':
    main()
