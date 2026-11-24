#!/usr/bin/env python
"""Write the test dataset test/adata_test.h5ad from a public Parse Biosciences dataset.

Source: "Performance of Evercode WT in Mouse Lymph Node Nuclei", Parse Biosciences (Evercode WT v2,
split-pipe v1.0.5), CC BY 4.0: https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-in-mouse-lymph-node-nuclei/

The download is a filtered matrix (called nuclei only), so the script:
  - takes a random subset of nuclei and splits them into two pseudo-plates (`plate`) and samples (`sample`);
  - keeps the genes detected in those nuclei, plus every mt- gene;
  - adds simulated empty droplets (10-100 counts drawn from the summed profile of all nuclei), so cell
    calling has background to separate. They are listed in `uns['test_data']`.

    python test/make_test_data.py --download-dir /tmp/parse_ln     # downloads ~310 MB once, checks sha256
"""
import argparse
import hashlib
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

HERE = Path(__file__).parent
SOURCE = 'https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-in-mouse-lymph-node-nuclei/'
BASE = 'https://cdn.parsebiosciences.com/lymph-node/2023/'
FILES = {
    'DGE.mtx': '40f9b25b6906d84c50e7f1267dd252ea8d1c960d6ec6792fe0dd8ab23b05401a',
    'all_genes.csv': '7200d60cd6a813c7b2c67399a38d4201953c7eef44c02499235bb120ef022db6',
    'cell_metadata.csv': 'd3e3c5d6d491ad046668e62f5b1a50024ec5dfb23abe795c8b76afac6e2c229b',
}
PLATES = {'plate1': ['LN_A', 'LN_B'], 'plate2': ['LN_C', 'LN_D']}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(download_dir):
    download_dir.mkdir(parents=True, exist_ok=True)
    for name, expected in FILES.items():
        path = download_dir / name
        if not path.exists():
            # download to a temporary name and rename only after the checksum matches,
            # so an interrupted download is retried on the next run
            tmp = path.with_suffix(path.suffix + '.part')
            print(f'downloading {BASE + name}')
            try:
                urllib.request.urlretrieve(BASE + name, tmp)  # noqa: S310 (fixed https URL)
                if sha256(tmp) != expected:
                    raise SystemExit(f'{BASE + name}: sha256 does not match the pinned source file')
                tmp.replace(path)
            finally:
                tmp.unlink(missing_ok=True)
        if sha256(path) != expected:
            raise SystemExit(f'{path}: sha256 does not match the pinned source file')
    return download_dir


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--download-dir', type=Path, required=True, help='where the Parse files are (or will be) saved')
    p.add_argument('--cells', type=int, default=300, help='nuclei per plate')
    p.add_argument('--empties', type=int, default=1500, help='simulated empty droplets per plate')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', type=Path, default=HERE / 'adata_test.h5ad')
    return p


def main():
    args = build_parser().parse_args()
    rng = np.random.default_rng(args.seed)
    src = fetch(args.download_dir)

    X = scipy.io.mmread(src / 'DGE.mtx').tocsr()                    # nuclei x genes
    obs = pd.read_csv(src / 'cell_metadata.csv', index_col='bc_wells')
    var = pd.read_csv(src / 'all_genes.csv').set_index('gene_name')
    var.index = pd.Index(ad.utils.make_index_unique(var.index.to_series().astype(str)))
    ambient = np.asarray(X.sum(0)).ravel()
    ambient = ambient / ambient.sum()

    n_plates = len(PLATES)
    pick = rng.choice(X.shape[0], size=args.cells * n_plates, replace=False)
    cells = X[pick]
    cell_obs = obs.iloc[pick].copy()

    keep = (np.asarray((cells > 0).sum(0)).ravel() > 0) | var.index.str.startswith('mt-')
    cells, ambient, var = cells[:, keep], ambient[keep], var.loc[keep]
    ambient = ambient / ambient.sum()

    depth = rng.integers(10, 100, size=args.empties * n_plates)
    empties = sp.csr_matrix(rng.poisson(ambient[None, :] * depth[:, None]))
    empty_obs = pd.DataFrame(index=[f'empty_{i:05d}' for i in range(len(depth))])
    empty_obs['tscp_count'] = np.asarray(empties.sum(1)).ravel()
    empty_obs['gene_count'] = np.asarray((empties > 0).sum(1)).ravel()
    empty_obs['mread_count'] = (empty_obs['tscp_count'] * rng.uniform(1.1, 1.6, size=len(depth))).astype(int)
    for bc in ('bc1', 'bc2', 'bc3'):
        well = rng.integers(1, 97, size=len(depth))
        empty_obs[f'{bc}_wind'] = well
        empty_obs[f'{bc}_well'] = [f'{"ABCDEFGH"[(w - 1) // 12]}{(w - 1) % 12 + 1}' for w in well]
    empty_obs['species'] = 'mm10'

    obs = pd.concat([cell_obs, empty_obs])
    plate_names = list(PLATES)
    plates = np.concatenate([np.repeat(plate_names, args.cells), np.repeat(plate_names, args.empties)])
    obs['plate'] = plates
    obs['sample'] = [PLATES[p][i % 2] for i, p in enumerate(plates)]
    obs['sublibrary'] = [f'{p}_s1' for p in plates]
    obs.index = [f'{bc}__{p}' for bc, p in zip(obs.index, plates)]
    obs = obs.drop(columns=[c for c in obs.columns if c not in
                            ['plate', 'sample', 'sublibrary', 'species', 'tscp_count', 'gene_count', 'mread_count',
                             'bc1_well', 'bc1_wind', 'bc2_well', 'bc2_wind', 'bc3_well', 'bc3_wind']])
    for col in ('plate', 'sample', 'species', 'sublibrary'):
        obs[col] = obs[col].astype('category')

    adata = ad.AnnData(X=sp.vstack([cells, empties]).tocsr().astype(np.float32), obs=obs, var=var[['gene_id']])
    adata.uns['test_data'] = {
        'source': SOURCE, 'license': 'CC BY 4.0', 'seed': args.seed,
        'nuclei': int(cells.shape[0]), 'simulated_empty_droplets': int(empties.shape[0]),
    }
    adata.write_h5ad(args.out, compression='gzip')
    print(f'wrote {args.out}: {cells.shape[0]} nuclei + {empties.shape[0]} empty droplets x {adata.n_vars} genes')


if __name__ == '__main__':
    main()
