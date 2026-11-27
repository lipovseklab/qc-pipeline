# qc-pipeline

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)

QC pipeline for mouse utricle single-cell data: cell calling (CellBender),
doublet detection (vaeda + scDblFinder), minimal filtering, dimensionality
reduction (TopOMetry), clustering (scSHC + clustree) and outlier removal; then
gene expression programs with cNMF on the annotated data.

Laid out following nf-core conventions (`main.nf`, `workflows/`, `modules/local/<name>/main.nf`, `conf/`,
`nextflow_schema.json`). `nf-core pipelines lint` has not been run yet.

## Running

Everything about a run goes in **one file**, `run.yaml`; copy one of the examples and edit it:

1. `cp examples/outliers.yaml run.yaml` (or `qc.yaml`, `all.yaml`, `cnmf.yaml`)
2. Edit `step`, `outdir`, `container_dir` (absolute path to the images, see below) and the dataset entry. Only
   settings that differ from the defaults are needed.
3. Run:
   ```bash
   nextflow run lipovseklab/qc-pipeline -r <version> -profile ucl_myriad -params-file run.yaml
   ```

`step` is one of:

| step | does | input `adata` |
|---|---|---|
| `qc` | cell calling (optional) → QC plots/stats → minimal filter → TopOMetry → clustering → doublets | raw / cell-called adata |
| `outliers` | MAD-based outlier removal per cluster → QC plots → filtering summary → TopOMetry → clustering | `qc/anndatas/adata_doublets.h5ad` from a qc run |
| `all` | `qc`, then `outliers` on its output | raw / cell-called adata |
| `cnmf` | cNMF per subset: subset → TopOMetry (HVGs) → prepare → factorize → combine + k selection plot; clustree | annotated adata (raw counts) |
| `consensus` | cNMF consensus for the chosen k + usage refit | same as `cnmf`, same outdir |
| `analysis` | GEP analysis: usages, top genes, enrichment, DGE, figures | same as `cnmf`, same outdir |

Outputs go to `<outdir>/qc/...`, `<outdir>/outliers/...` and `<outdir>/cnmf/<subset>/...`. The cNMF steps are
separate runs because the k and the clusters are chosen from plots in between: see **[docs/cnmf.md](docs/cnmf.md)**
for their settings (a `cnmf:` block under the dataset) and outputs.

The run file is checked before anything starts: misspelt settings, wrong types (`skip_scSHC: "yes"`),
missing required settings and missing input files all stop the run with a message saying what to fix.

### Provenance: every result is tied to a commit

- The pipeline **stops if its code has uncommitted changes**, so results can always be traced back to the exact code.
  Running from GitHub with `-r <tag>` always passes. For a test run from an edited local copy, add `--allow_dirty true`.
- Each run writes `<outdir>/pipeline_info/provenance.json`: commit, tag (`git describe`), the edited files
  (if `--allow_dirty`), run name and session ID, command line, profile, Nextflow version, the dataset settings
  after defaults, all parameters, the final status, and the sha256 of each container image in `container_dir`
  (cached next to the image as `<image>.sif.sha256`).

## Dataset settings

One entry under `datasets:` (one dataset per run for now).

| setting | needed for | default | meaning |
|---|---|---|---|
| `name` | all | — | short name, used in output file names (`alldata`) |
| `adata` | all | — | input AnnData (see table above) |
| `celltype_csv` | qc, outliers, all | — | marker table for simple cell typing (format below) |
| `batch` | all | `plate` | obs column with the batch |
| `sample` | all | `sample` | obs column used to colour QC plots |
| `skip_scSHC` | all | `false` | skip scSHC clustering (needs ~128 GB) |
| `skip_cellbender` | qc | `true` | `adata` is already cell-called (see below) |
| `cellbender_csv` | qc | — | per-batch CellBender settings: `batch, total_droplets_included, expected_cells` (example: `examples/cellbender.csv`) |
| `clusters` | outliers, all | — | obs column with the clusters MADs are computed in, e.g. `topo_clusters_ms_res0.4` |
| `doublet_column` | outliers, all | `scDblFinder_calls` | obs column with `doublet`/`singlet` calls |
| `qc_stats` | outliers | — | the qc run's `qc/reports/stats` folder, so the filtering summary covers every stage |
| `cnmf` | cnmf, consensus, analysis | — | cNMF settings and subsets ([docs/cnmf.md](docs/cnmf.md)) |

Run-wide settings: `outdir` (default `results`) and `container_dir` (folder with the `.sif` images; default
`<pipeline>/containers`; `~/Scratch/containers` with `-profile ucl_myriad`). The images are not in git: build them once
with `bin/helpers/setup-containers.sh` (in a clone of this repo; `--dest <folder>` puts them elsewhere, e.g.
`--dest ~/Scratch/containers` on Myriad) and set `container_dir` to the **absolute** path of that folder in
`run.yaml` unless the profile already points there. `nextflow run lipovseklab/qc-pipeline` runs its own cached copy of the code, which has no images.

**Cell calling.** The default `skip_cellbender: true` expects an already **cell-called** `adata` (nuclei with
`cell_probability < 0.5` removed), e.g. from a separate CellBender run; no empty-droplet removal happens in this
mode. With `skip_cellbender: false`, `adata` must be the fully unfiltered matrix, and the machine needs an NVIDIA
GPU: `run_cellbender` runs with `--cuda` and the container gets `--nv`. vaeda runs on the CPU (as in the paper run;
on the paper data a GPU run ran out of GPU memory).

**`celltype_csv`** has **no header row**: column 1 is a cell-type label, the remaining columns are marker gene
symbols (rows can differ in length; see `test/celltype_test.csv`). Markers not in `adata.var_names` are ignored,
so check gene symbols match your naming.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 25.10 (the version CI tests is 25.10.4; nf-schema 2.5.1 and workflow outputs need it)
- [Apptainer/Singularity](https://apptainer.org/) — all processes run in prebuilt containers.

## Profiles

Profiles hold machine settings only; dataset settings always go in `run.yaml`.

- `test` — small public dataset (`test/adata_test.h5ad`: 600 Parse mouse lymph node nuclei plus simulated empty droplets, CC BY 4.0; see `test/README.md`), step `all`; used by `nf-test`.
- `ucl_myriad` — UCL Myriad on its new operating system and Slurm scheduler (log in to `myriad-ng.rc.ucl.ac.uk`;
  the old SGE half is being retired). Images in `~/Scratch/containers` unless `run.yaml` sets `container_dir`.
  Nextflow needs Java 17 or later and the system Java is 8, so load it first, then run Nextflow from a login node
  inside `tmux` or `screen` so it survives a dropped connection:
  ```bash
  module load openjdk/21.0.3_9/gcc-12.3.0-kx7sucd
  tmux new -s qc
  nextflow run lipovseklab/qc-pipeline -r <version> -profile ucl_myriad -params-file run.yaml
  ```
  CellBender jobs go to the `gpu` partition with `--gres=gpu:1`; the other steps get memory and time sized from the
  paper runs (`conf/ucl_myriad.config`).

## Testing

```bash
nf-test test                                   # stub runs: qc (with and without CellBender), outliers, all,
                                               # cnmf, consensus, analysis, and the run.yaml checks
nextflow lint . -exclude nf-test.config        # syntax check (needs Nextflow >= 25.04; CI uses 25.10.4)
```

Stub runs need no containers or data. Both checks run on every pull request (`.github/workflows/ci.yml`).

## License

Released under the standard MIT license.
