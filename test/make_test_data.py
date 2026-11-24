#!/usr/bin/env python
"""Write the synthetic test dataset test/adata_test.h5ad.

The counts are random (no real data): a few made-up cell types, each with its own gene profile and higher
counts for its markers from test/celltype_test.csv, plus low-count empty droplets so cell calling has
something to separate. obs columns match a split-pool (Parse) unfiltered matrix.

    python test/make_test_data.py            # writes test/adata_test.h5ad
"""
import argparse
import csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).parent
PLATES = {'plate1': ['E14', 'P0'], 'plate2': ['P7', 'P21']}


def marker_genes(path):
    with open(path) as fh:
        rows = list(csv.reader(fh))[1:]
    return {row[0]: [g for g in row[1:] if g] for row in rows}


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--cells', type=int, default=300, help='real nuclei per plate')
    p.add_argument('--empties', type=int, default=1500, help='empty droplets per plate')
    p.add_argument('--genes', type=int, default=2000, help='total genes')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', default=HERE / 'adata_test.h5ad')
    return p


def main():
    args = build_parser().parse_args()
    rng = np.random.default_rng(args.seed)

    markers = marker_genes(HERE / 'celltype_test.csv')
    known = sorted({g for genes in markers.values() for g in genes})
    mito = [f'mt-{g}' for g in ['Nd1', 'Nd2', 'Co1', 'Co2', 'Atp8', 'Atp6', 'Co3', 'Nd3', 'Nd4l', 'Nd4', 'Nd5', 'Nd6', 'Cytb']]
    filler = [f'Gene{i:04d}' for i in range(args.genes - len(known) - len(mito))]
    genes = known + mito + filler
    gene_index = {g: i for i, g in enumerate(genes)}

    # one expression profile per cell type, markers boosted
    types = list(markers)
    base = rng.gamma(0.3, 1.0, size=len(genes))
    profiles = []
    for t in types:
        prof = base * rng.lognormal(0, 0.5, size=len(genes))
        for g in markers[t]:
            prof[gene_index[g]] *= 50
        profiles.append(prof / prof.sum())
    ambient = base / base.sum()

    blocks, obs = [], []
    for plate, samples in PLATES.items():
        for kind, n in (('cell', args.cells), ('empty', args.empties)):
            if kind == 'cell':
                depth = rng.lognormal(np.log(2500), 0.5, size=n).astype(int) + 500
                probs = np.array([profiles[i] for i in rng.integers(len(types), size=n)])
            else:
                depth = rng.integers(10, 100, size=n)
                probs = np.tile(ambient, (n, 1))
            counts = rng.poisson(probs * depth[:, None]).astype(np.float32)
            blocks.append(sp.csr_matrix(counts))
            for i in range(n):
                obs.append({'plate': plate, 'sample': samples[i % len(samples)], 'sublibrary': f'{plate}_s1'})
    X = sp.vstack(blocks).tocsr()

    obs = pd.DataFrame(obs)
    obs.index = [f'bc{i:06d}__{p}' for i, p in enumerate(obs['plate'])]
    obs['species'] = 'M'
    obs['tscp_count'] = np.asarray(X.sum(1)).ravel().astype(int)
    obs['gene_count'] = np.asarray((X > 0).sum(1)).ravel().astype(int)
    obs['mread_count'] = (obs['tscp_count'] * rng.uniform(1.1, 1.6, size=len(obs))).astype(int)
    for bc in ('bc1', 'bc2', 'bc3'):
        well = rng.integers(1, 97, size=len(obs))
        obs[f'{bc}_well'] = [f'{"ABCDEFGH"[(w - 1) // 12]}{(w - 1) % 12 + 1}' for w in well]
        obs[f'{bc}_wind'] = well
    for col in ('plate', 'sample', 'species', 'sublibrary'):
        obs[col] = obs[col].astype('category')

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.write_h5ad(args.out, compression='gzip')
    print(f'wrote {args.out}: {adata.n_obs} barcodes x {adata.n_vars} genes')


if __name__ == '__main__':
    main()
