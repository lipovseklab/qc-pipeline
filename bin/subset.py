#!/usr/bin/env python3
# import packages
import scanpy as sc
import anndata as ad
import argparse
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))   # bin/ modules, wherever this runs from
from plain_strings import plain_strings, report

def build_parser():
	
	p = argparse.ArgumentParser(
        	description=__doc__,
        	formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	p.add_argument(
		"--counts", 
		metavar="FILE",
                help="AnnData file (.h5ad)"
	)
	p.add_argument(
		"--metadata", 
                help="Metadata column in obs e.g. celltype"
	)
	p.add_argument(
		"--subset", 
        	help="Subset of metadata to keep e.g. HC"
		)
	p.add_argument(
		"--label",
		default=None,
		help="Name for the output file adata_{metadata}_{label}.h5ad (default: --subset)"
	)
	p.add_argument(
		'--output-dir',
                help="top-level directory for saving outputs"
	)
    
	return p

def main():
	
	# setup arguments
	args = build_parser().parse_args()
	
	# setup outpath
	outpath = Path(args.output_dir)
	outpath.mkdir(parents=True, exist_ok=True)
	# testing
	print(outpath)
	
	# read adata
	adata = ad.read_h5ad(args.counts)
	report(args.counts, plain_strings(adata))

	# subset by obs; the value arrives as text, so compare as a number when the column is numeric
	column = adata.obs[args.metadata]
	value = args.subset
	dtype = column.cat.categories.dtype if isinstance(column.dtype, pd.CategoricalDtype) else column.dtype
	if pd.api.types.is_numeric_dtype(dtype):
		value = pd.to_numeric(value)
	# keep every category of the parent (e.g. all ages, not only those in this subset), so colours made later cover
	# the full list and a group has the same colour in every subset; stored colours (uns['<col>_colors']) stay aligned
	# plain text columns become categories too (e.g. an age column stored as strings), so they keep every value
	for col in adata.obs.columns:
		if pd.api.types.infer_dtype(adata.obs[col], skipna=True) == 'string':
			adata.obs[col] = adata.obs[col].astype('category')
	with ad.settings.override(remove_unused_categories=False):
		adata = adata[column == value].copy()

	# write adata
	adata.write(f'{outpath}/adata_{args.metadata}_{args.label or args.subset}.h5ad')

if __name__ == "__main__":
    main()

