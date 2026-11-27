#!/usr/bin/env python3
"""Write an AnnData with raw counts in X, for `cnmf prepare` (which reads X as counts).

Takes .raw if present, else layers['counts'], else X when it already holds whole numbers; stops otherwise.
Used for whole-dataset cNMF runs given an HVG list (other runs get their counts from topometry.py).

    cnmf_counts.py --adata in.h5ad --out counts.h5ad
"""
import argparse

import anndata as ad
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))   # bin/ modules, wherever this runs from
from plain_strings import plain_strings, report
import numpy as np
import scipy.sparse as sp


def whole_numbers(X):
    values = X.data[:1000] if sp.issparse(X) else np.asarray(X).ravel()[:1000]
    return np.allclose(values, np.round(values), atol=1e-6)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--adata', required=True)
    p.add_argument('--out', required=True)
    return p


def main():
    args = build_parser().parse_args()
    adata = ad.read_h5ad(args.adata)
    report(args.adata, plain_strings(adata))
    if adata.raw is not None:
        counts, source = adata.raw.to_adata(), '.raw'
    elif 'counts' in adata.layers:
        counts, source = ad.AnnData(X=adata.layers['counts'], obs=adata.obs, var=adata.var), "layers['counts']"
    else:
        counts, source = adata, 'X'
    if not whole_numbers(counts.X):
        raise SystemExit(f'{args.adata}: no raw counts found (.raw, layers["counts"] or X with whole numbers)')
    print(f'{args.adata}: counts from {source}, {counts.n_obs} cells x {counts.n_vars} genes')
    counts.write_h5ad(args.out)


if __name__ == '__main__':
    main()
