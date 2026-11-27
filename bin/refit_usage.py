#!/usr/bin/env python3
"""Repeat the final usage refit of `cnmf consensus` with the genes in matching order, and check cNMF's own usages.

cNMF 1.7.1 finishes `consensus` (default refit_usage=True) by regressing the std-scaled TPM of the overdispersed genes on
the consensus spectra (TPM units, divided by the same std):

    spectra_tpm_rf = spectra_tpm.loc[:, hvgs].div(tpm_stats.loc[hvgs, '__std'], axis=1)

`spectra_tpm` is labelled with the gene names of `tpm.h5ad`, `tpm_stats` with plain strings. When the names in
`tpm.h5ad` are pandas' nullable `string` type, pandas treats the two labels as different, aligns them by sorting the
genes alphabetically, and the next step pairs the spectra and the data by position on different genes: the usages in
`*.usages.k_*.consensus.txt` are then noise. With plain names both agree.

This runs the same step with cNMF's own `refit_usage` and NMF settings on plain arrays, writes the result next to
the cNMF files as `<name>.usages.k_<k>.dt_<threshold>.consensus.refit.txt`, and stops (exit 1, after writing) if
cNMF's own usages differ from it (`--min_r`).

    refit_usage.py --output_dir results/ctype_MEL --name cnmf_run --k 40 --threshold 0.2
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from cnmf import cNMF


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output_dir', required=True, help='folder holding the cNMF run folder (as given to `cnmf --output-dir`)')
    p.add_argument('--name', required=True, help='cNMF run name')
    p.add_argument('--k', type=int, required=True, help='number of components given to `cnmf consensus`')
    p.add_argument('--threshold', type=float, required=True, help='local density threshold given to `cnmf consensus`')
    p.add_argument('--top_genes', type=int, default=10, help='top genes per GEP used for the sanity check')
    p.add_argument('--min_r', type=float, default=0.99,
                   help="stop if any GEP's refit usage correlates below this with cNMF's own usage (0 or less: only report)")
    return p.parse_args()


def main():
    args = parse_args()
    run_dir = Path(args.output_dir) / args.name
    dt = str(float(args.threshold)).replace('.', '_')          # how cNMF writes the threshold into file names
    c = cNMF(output_dir=args.output_dir, name=args.name)

    spectra_file = run_dir / f'{args.name}.gene_spectra_tpm.k_{args.k}.dt_{dt}.txt'
    score_file = run_dir / f'{args.name}.gene_spectra_score.k_{args.k}.dt_{dt}.txt'
    for f in (spectra_file, score_file, Path(c.paths['tpm']), Path(c.paths['tpm_stats']),
              Path(c.paths['nmf_genes_list']), Path(c.paths['nmf_run_parameters'])):
        if not f.exists():
            raise SystemExit(f'{f} not found. Has `cnmf consensus` been run for k={args.k}, threshold={args.threshold}?')

    genes = [g for g in open(c.paths['nmf_genes_list']).read().split('\n') if g]
    spectra = pd.read_csv(spectra_file, sep='\t', index_col=0)           # GEPs x all genes, TPM units
    with np.load(c.paths['tpm_stats'], allow_pickle=True) as f:
        std = pd.DataFrame(**f).loc[genes, '__std'].to_numpy()

    tpm = ad.read_h5ad(c.paths['tpm'])
    data = tpm[:, genes].copy()
    if sp.issparse(data.X):
        sc.pp.scale(data, zero_center=False)
    else:
        data.X = data.X / data.X.std(axis=0, ddof=1)
    assert list(data.var_names) == genes

    # both matrices in `genes` order: plain arrays, so nothing can re-sort the genes by label
    h = pd.DataFrame((spectra.loc[:, genes].to_numpy() / std).astype(data.X.dtype), index=spectra.index, columns=genes)
    usage = pd.DataFrame(np.asarray(c.refit_usage(data.X, h)), index=tpm.obs_names, columns=spectra.index)

    if usage.index.duplicated().any():
        dups = usage.index[usage.index.duplicated()].unique()
        raise SystemExit(f'{len(dups)} cell IDs occur more than once (e.g. {list(dups[:3])}); cNMF needs unique cell '
                         f'IDs. Make them unique in the input AnnData (adata.obs_names_make_unique()) and rerun.')
    out = run_dir / f'{args.name}.usages.k_{args.k}.dt_{dt}.consensus.refit.txt'
    usage.to_csv(out, sep='\t')
    print(f'wrote {out} ({usage.shape[0]:,} cells x {usage.shape[1]} GEPs)')

    # check: cNMF's own usages should be the same as this refit. They differ when the run folder was prepared from
    # gene names in pandas' `string` type (cNMF then pairs genes in alphabetical order; docs/cnmf.md)
    own_file = run_dir / f'{args.name}.usages.k_{args.k}.dt_{dt}.consensus.txt'
    # cell IDs as text: numeric IDs, leading zeros and IDs such as 'NA' stay as written
    own = pd.read_csv(own_file, sep='\t', index_col=0, dtype=str, keep_default_na=False)
    own.columns = own.columns.astype(str)
    own = own.loc[usage.index.astype(str), [str(c) for c in usage.columns]].astype(float).to_numpy()
    same = pd.Series({gep: np.corrcoef(usage[gep].to_numpy(), own[:, j])[0, 1] for j, gep in enumerate(usage.columns)})
    print(f'r(refit usage, cNMF usage) per GEP: median {same.median():.3f}, min {same.min():.3f}')
    if args.min_r > 0 and not (same > args.min_r).all():
        raise SystemExit(
            f'cNMF\'s own usages ({own_file.name}) differ from the refit (min r {same.min():.2f} < {args.min_r}). '
            f'The run folder was probably prepared before names were stored as plain strings: its *.consensus.txt '
            f'usages are wrong and {out.name} is right. Rerun the cnmf step (prepare) to rebuild the folder, or rerun '
            f'with --min_r 0 to accept the refit.')

    # sanity check: usage should follow the expression of the GEP's top genes
    scores = pd.read_csv(score_file, sep='\t', index_col=0)              # GEPs x genes
    top = {gep: scores.loc[gep].sort_values(ascending=False).index[:args.top_genes] for gep in scores.index}
    need = sorted({g for v in top.values() for g in v})
    expr = tpm[:, need].X
    expr = pd.DataFrame(np.log1p(expr.toarray() if sp.issparse(expr) else expr), index=tpm.obs_names, columns=need)
    r = pd.Series({gep: np.corrcoef(expr[list(top[gep])].mean(axis=1), usage[gep])[0, 1] for gep in scores.index})
    print(f'r(usage, mean expression of the top {args.top_genes} genes) per GEP: median {r.median():.2f}, '
          f'range {r.min():.2f} to {r.max():.2f}')
    if r.median() < 0.3:
        print('WARNING: the usages hardly follow the top genes; check the cNMF outputs', file=sys.stderr)


if __name__ == '__main__':
    main()
