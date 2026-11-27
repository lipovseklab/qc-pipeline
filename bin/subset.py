#!/usr/bin/env python3
# import packages
import scanpy as sc
import anndata as ad
import argparse
import pandas as pd
from pathlib import Path

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

	# subset by obs; the value arrives as text, so compare as a number when the column is numeric
	column = adata.obs[args.metadata]
	value = args.subset
	dtype = column.cat.categories.dtype if isinstance(column.dtype, pd.CategoricalDtype) else column.dtype
	if pd.api.types.is_numeric_dtype(dtype):
		value = pd.to_numeric(value)
	adata = adata[column == value]

	# write adata
	adata.write(f'{outpath}/adata_{args.metadata}_{args.label or args.subset}.h5ad')

if __name__ == "__main__":
    main()

