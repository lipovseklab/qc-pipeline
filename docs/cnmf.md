# cNMF: gene expression programs

The cNMF steps run [cNMF](https://github.com/dylkot/cNMF) on an annotated dataset, on the whole dataset or on
subsets of it (e.g. one cell type each), and then analyse the gene expression programs (GEPs). They are three
runs on the same `outdir`, because two choices have to be made by looking at plots in between:

| step | does | you then choose |
|---|---|---|
| `cnmf` | per subset: subset → TopOMetry (HVGs, clusterings) → `cnmf prepare` → `factorize` → `combine` + k selection plot; cluster comparison + clustree plots | `selected_k` and `threshold` from the k selection plot |
| `consensus` | `cnmf consensus` for the chosen k + usage refit (below) | `clusters` from the clustree plots |
| `analysis` | GEP usages into the AnnData, top genes, GO BP + CollecTRI enrichment, DGE by cluster, figures | |

Copy `examples/cnmf.yaml`, edit it, and run each step with the same file, changing only `step` and adding the
choices (`selected_k`, `threshold`, `clusters`) to the subsets. The file then records every choice of the run.

```bash
nextflow run lipovseklab/qc-pipeline -r <version> -profile ucl_myriad -params-file run.yaml
```

Run in `tmux` or `screen` on a login node: Nextflow submits the jobs and has to stay alive meanwhile.

## Settings

Under the dataset in `run.yaml`:

```yaml
datasets:
  - name: alldata
    adata: /path/to/adata_ctype.h5ad
    cnmf:
      column: ctype            # settings here apply to every subset ...
      k: [20, 30, 40]
      subsets:
        - value: HC            # ... and can be overridden per subset
          k: [10, 15, 20]
          selected_k: 15
        - value: SC
```

- Each entry under `subsets` is one cNMF run on the cells with `obs[column] == value`, in the folder
  `<outdir>/cnmf/<column>_<value>` (characters other than letters, digits, `.`, `_` and `-` in the value become
  `_`, so `Hair cell` gives `ctype_Hair_cell`).
- An entry without `value` (`- {}`), or no `subsets` at all, runs on the **whole dataset**, in
  `<outdir>/cnmf/<dataset name>`. TopOMetry selects its HVGs, unless the entry gives `hvgs:` (a gene list, one per
  line), which skips TopOMetry and starts from `cnmf prepare` on `adata`.
- The `adata` needs raw counts in `.raw` (or `layers['counts']`); a subset with 0 cells stops the run.

| setting | step | default | meaning |
|---|---|---|---|
| `column` | cnmf | — | obs column the subsets are taken from (needed with `value`) |
| `value` | cnmf | — | per subset: the value of `column` to keep |
| `k` | cnmf | `[20, 30, 40]` | numbers of components to factorize |
| `n_iters` | cnmf | `100` | iterations (replicates) per k |
| `stride` | cnmf | `10` | factorize iterations per job |
| `seed` | cnmf | `42` | `cnmf prepare` seed; a run folder made with another seed stops the run (no silent mixing) |
| `hvgs` | cnmf | — | whole dataset only: gene list to use instead of TopOMetry's HVGs |
| `run_name` | all | see below | cNMF run name |
| `cluster_color` | cnmf | `ctype` | obs column shown next to the clusterings in the comparison plot (skipped if missing) |
| `selected_k` | consensus, analysis | — | the k chosen from the k selection plot; subsets without it are skipped |
| `threshold` | consensus, analysis | `2` | local density threshold (`2` = no filtering) |
| `clusters` | analysis | `clusters` | obs column with the clusters, chosen from the clustree plots, e.g. `topo_clusters_ms_res0.6` |
| `projection` | cnmf, analysis | `projection` | obsm key of the embedding to plot on |
| `age` | analysis | `age_pretty` | obs column with the age labels |
| `summary_cols` | analysis | age, clusters, `ctype_detailed` | columns for the projection summary (missing ones are skipped) |
| `usage_cutoff` | analysis | `0.1` | usage above which a cell counts as using a GEP |
| `n_top_genes` | analysis | `100` | top genes per GEP |
| `geps_per_page` | analysis | `6` | GEPs per page in the multi-GEP figures (`0` = one page) |
| `formats` | analysis | `png,pdf` | figure formats |
| `analysis_dir` | analysis | `analysis` | output folder name; use another to keep a second analysis of the same subset |

**Run name.** The `cnmf` step names each run `cnmf_<tag>_<YYYYMMDD>` and records it in
`<outdir>/cnmf/<tag>/pipeline_info/run_name.txt`. Later runs into the same outdir (a resume on another day,
`consensus`, `analysis`) read it from there; give `run_name` to override it.

**Resuming.** `-resume` continues the latest run in the launch folder; to continue a specific run, give its session
ID (`nextflow log`), not its run name. Don't resume while one of its jobs is still running (`qstat`).

## Usage refit

cNMF 1.7.1 ends `consensus` with a refit of the usages, but pairs the spectra (genes sorted alphabetically) with the
data (genes in HVG order) by position, so its `*.usages.k_*.consensus.txt` is wrong: the usages do not follow each
GEP's top genes (spectra, gene scores and top genes are fine). The `consensus` step therefore runs
`bin/refit_usage.py` afterwards, which repeats the refit with cNMF's own function and settings and both matrices in
the same gene order, and writes `<run_name>.usages.k_<k>.dt_<threshold>.consensus.refit.txt` next to cNMF's files
(which are left as they are). It prints how well each GEP's usage follows the mean expression of its top 10 genes
(about 0.5–0.9 when right, about 0 for cNMF's own file). The analysis reads the refitted file.

## Analysis

One-off setup on a machine with internet access (compute nodes may have none), in a clone of this repo:

```bash
apptainer exec containers/cnmf-analysis.sif python bin/fetch_resources.py   # writes assets/resources/ (GMT + CollecTRI)
```

`nextflow run lipovseklab/qc-pipeline` runs its own cached copy of the code, which does not have these files, so
give their **absolute** paths in `run.yaml` (top level, next to `container_dir`):

```yaml
gmt: /path/to/qc-pipeline/assets/resources/MOUSE_GO_bp_no_GO_iea_symbol.gmt
collectri: /path/to/qc-pipeline/assets/resources/collectri_mouse.parquet
```

Outputs go to `<outdir>/cnmf/<tag>/analysis/`:

```
csvs/       top100_genes.csv, gsea_gobp.csv, gsea_collectri.csv, dge_by_cluster.csv
excel/      top100_genes.xlsx, gsea_gobp.xlsx, gsea_collectri.xlsx   (gsea: one sheet per GEP)
figures/    gep_usage_pct_by_cluster, gep_usage_pct_by_age, gep_usage_heatmap_cluster_age, projection_summary,
            gep_megaplot_GEP01-06, gep_usage_topgenes_GEP01-06, gep_gobp_GEP01-06, gep_collectri_GEP01-06, ...
            elements/   every panel on its own
anndatas/   adata_cnmf_<tag>.h5ad                                    (usages in obs, gene scores in varm)
analysis_call.txt, analysis_call.json                                 (how it was run)
```

- `analysis_call.txt` / `.json` record the `nextflow run` command, run name and session ID, the `gep_analysis.py`
  command with every option written out, the input files (paths, sizes, dates) and the package versions.
- Gene panels and the DGE use log-normalised expression of all genes (genes in ≥ 3 cells, `normalize_total` to
  1e4, `log1p`). `dge_by_cluster.csv`: each cluster against the rest (Wilcoxon, `tie_correct=True`); clusters with
  fewer than 10 cells are left out.
- Figures are PNG and PDF (points rasterised, text editable). Age uses seaborn's `husl` palette; each cluster column
  gets one palette used by every figure, stored in the output AnnData (`uns['<column>_colors']`).
- `notebooks/gep_analysis_template.ipynb` lays out the same steps one by one with the functions in `bin/gep_lib`.

## Outputs

```
<outdir>/cnmf/<tag>/
    <run_name>/                 cNMF run folder (k selection plot, spectra, usages, consensus files)
    anndatas/                   adata_<tag>.h5ad (subset), adata_topometry_<tag>.h5ad
    hvgs_<tag>.csv              the HVGs given to cNMF
    topometry/                  TopOMetry plots and object
    clustering/                 clustree_<series>.png, cluster_comparison_on_embedding_<series>.png
    analysis/                   see above
    pipeline_info/run_name.txt
<outdir>/pipeline_info/         provenance.json, execution trace/report/timeline of each run
```

The clustree and cluster comparison plots are made for every clustering series in the TopOMetry AnnData — `pca`
(`pca_leiden_res*`), `topo` (`topo_clusters_res*`) and `topo_ms` (`topo_clusters_ms_res*`) — so the latent spaces can
be compared before choosing `clusters`.

## Containers

`cnmf_env.sif` (cNMF 1.7.1 from BioContainers), `cnmf-analysis.sif`, `clustree.sif`, and `utricle-qc.sif` for
TopOMetry, all in `container_dir`. Build them into that folder from a clone of this repo, e.g. on Myriad (where the
`ucl_myriad` profile sets `container_dir` to `~/Scratch/containers`):

```bash
bin/helpers/setup-containers.sh --dest ~/Scratch/containers cnmf_env cnmf-analysis clustree utricle-qc
```

Without `--dest` the images go to the clone's `containers/`; then set `container_dir` to that folder in `run.yaml`.

## Resources

Set per label in `conf/base.config` (`cnmf_small`, `cnmf_topometry`, `cnmf_factorize`, `cnmf_combine`,
`cnmf_consensus`, `cnmf_analysis`), sized from the 2026 Myriad runs on the paper subsets (7–13k cells). Failed jobs
killed for memory or time are retried once with twice the resources. A much larger subset may need more.
