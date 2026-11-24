# Test data

`adata_test.h5ad` is made by `make_test_data.py` from a public dataset:

> Parse Biosciences, *Performance of Evercode™ WT in Mouse Lymph Node Nuclei* (Evercode WT v2, split-pipe v1.0.5),
> <https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-in-mouse-lymph-node-nuclei/>,
> licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Changes from the source: 600 random nuclei, split into two pseudo-plates (`plate1`, `plate2`) and four pseudo-samples;
genes not detected in those nuclei are dropped (mt- genes kept); 3,000 simulated empty droplets (10–100 counts drawn
from the summed profile of all nuclei) are added, because the download contains called nuclei only. Details are in
`adata.uns['test_data']`.

`celltype_test.csv` lists lymph node marker genes for the cell-type step.

Rebuild with `python test/make_test_data.py --download-dir <folder>` (downloads ~310 MB once and checks the
pinned sha256 of each file).
