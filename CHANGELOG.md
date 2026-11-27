# utricleqc: Changelog

## [Unreleased]

### `Fixed`

- cNMF usages: when the input stored gene names in pandas' nullable `string` type, cNMF 1.7.1's `consensus` paired
  genes in alphabetical order in its final usage refit, so `*.consensus.txt` usages were noise (the refit file was
  right). The subset, TopOMetry, counts and analysis scripts now store names and text columns (including `.raw`) as
  plain strings when they read their input (`bin/plain_strings.py`), and the refit step stops if cNMF's own usages
  differ from it (`refit_usage.py --min_r`, default 0.99).
- Test data: `test/adata_test.h5ad` is now built by `test/make_test_data.py` from a public Parse Biosciences mouse lymph
  node nuclei dataset (CC BY 4.0) plus simulated empty droplets; `test/celltype_test.csv` lists lymph node markers.
  `test/generate_test_from_raw.py` is removed.
- run.yaml checks: blank required settings count as missing; dataset `name` may only use letters, digits, `_`, `.`
  and `-`; text settings may not contain spaces or line breaks; a misspelt top-level setting (e.g. `outdri`) now stops the run.
- Provenance guard: a copy of the pipeline that isn't tracked by the git repository around it, or a failed
  `git status`, is no longer taken for a clean commit. An unreadable `.sha256` cache or an image-hashing error no
  longer stops the final `provenance.json` from being written.
- `bin/helpers/setup-containers.sh` keeps images that already exist (`--force` rebuilds), so it can be rerun after a
  failure, and can build single images (`setup-containers.sh utricle-qc`).
- README: minimum Nextflow is 25.10 (nf-schema 2.5.1 needs >= 25.04); `skip_cellbender: false` needs a GPU;
  `container_dir` must be an absolute path when running `lipovseklab/qc-pipeline` from GitHub.
- cNMF: a rerun with another `--seed` in an existing run folder no longer mixes old and new spectra: PREP records the
  seed in `<run>/cnmf_seed.txt` and stops on a mismatch (or on factorize results with no seed record).
- cNMF: TopOMetry PCA asks for at most `min(100, cells - 1, genes - 1)` components, so small subsets no longer fail.
- cNMF: run and sample names are limited to letters, digits, `_`, `.` and `-` (they are used in shell commands).
- cNMF: the Myriad profile uses `apptainer`, like the QC pipeline, so the two configs no longer enable two engines.
- cNMF README: `apptainer pull <file> <URI>` argument order and the correct image tag (`cnmf:1.7.1--pyhdfd78af_0`).
- `containers/utricle-qc.def` had a typo (`reBootstrap`) that stopped it building.

### `Changed`

- `ucl_myriad` profile: Myriad's new operating system and Slurm scheduler (`myriad-ng`) instead of SGE. QC and
  outlier steps get memory and time sized from the paper runs (Slurm's default is 1 GB per core), CellBender goes to
  the `gpu` partition with `--gres=gpu:1`, and steps without a time limit get 24 h. README: load Java 17+ on Myriad.
- cNMF is part of the main pipeline: `step: cnmf`, `consensus` and `analysis`, driven by the same `run.yaml` as QC
  (a `cnmf:` block under the dataset, with one entry per subset; `examples/cnmf.yaml`, `docs/cnmf.md`). This replaces
  `analysis/cnmf/` with its own `main.nf`, configs, samplesheet and `--sample`/`--selected_k` flags: the chosen k,
  threshold and clusters are now written per subset in `run.yaml`. Outputs go to `<outdir>/cnmf/<tag>/`. The images
  come from `container_dir` like the QC images (`setup-containers.sh cnmf_env cnmf-analysis clustree`).
- cNMF: `--step refit_usage` is gone; the consensus step always refits (rerun it to refit). Task logs are no longer
  copied to `<tag>/logs/` and the per-sample `trace.txt` is not kept: Nextflow's execution trace and report in
  `pipeline_info/` cover every run. A run name recorded in the outdir is reused, so `--run_date` is not needed to
  resume on another day.
- One `bin/topometry.py` for QC and cNMF: the cNMF copy's extra outputs (gene filter before HVG selection, HVG list,
  counts) are behind `--cnmf_outputs`. The PCA size guard for small subsets now applies to QC too (no change for
  datasets with more than 100 cells and HVGs). The cNMF clustree script is `bin/clustree_series.R`.
- vaeda runs on the CPU (no `--nv`): the paper run did (its TensorFlow found no CUDA libraries, and its calls match a
  CPU rerun exactly), and on the paper data a GPU run on a 5 GB GPU ran out of memory after 40 min.
- `ucl_myriad` profile: `container_dir` defaults to `~/Scratch/containers`; `queueSize` 800 for the cNMF factorize jobs.
- A stray word on the command line (e.g. `threshold 0.5` without `--`) stops the run, for every step.
- CI lints the whole pipeline.
- `bin/helpers/setup-containers.sh` passes `APPTAINER_BUILD_OPTS` to `apptainer build` (on WSL:
  `APPTAINER_BUILD_OPTS=--ignore-subuid` when a build stops with "failed to mount binfmt_misc filesystem").
- `bin/helpers/setup-containers.sh --dest <folder>` puts the images where `container_dir` points (e.g.
  `~/Scratch/containers`, the `ucl_myriad` default) instead of the clone's `containers/`.
- cNMF: a numeric subset value matches a numeric (or numeric categorical) obs column; `summary_cols` reach the
  analysis as one comma-separated list; a whole-dataset run with `hvgs` gives cNMF the raw counts (`.raw` or
  `layers['counts']`, via `bin/cnmf_counts.py`) instead of `X`; subsets without `selected_k` need no earlier run for
  consensus/analysis; `topometry.py --min_cells` sets both gene filters; the factorize and combine retry exceptions
  are no longer overridden by the common cNMF retry rule.

- Container recipes: removed the WSL `ssl_verify false` workaround; pinned `cellbender==0.3.0`, `scikit-misc==0.5.2`
  and `hnswlib` to the versions in the image used for the paper runs; dropped machine-specific `prefix:` lines.
- `bin/helpers/setup-containers.sh` works from any folder and stops on the first error.
- `params/cellbender_params.csv` moved to `examples/cellbender.csv`.

### `Removed`

- `analysis/cnmf/scripts/` (the older qsub route; the Nextflow pipeline replaces it) and the duplicate
  `analysis/cnmf/envs/utricle-qc.{def,yml}` (the QC pipeline's `containers/utricle-qc.def` is the one recipe).
- `bin/helpers/myriad.sh` (called a script that no longer exists), `bin/helpers/nextflow-rsync.sh` (replaced by the
  run archiving in the project repo) and the unused `params/join_params.csv`.

### `Added`

- Provenance guard: the pipeline stops on uncommitted code unless `--allow_dirty true`, and every run writes
  `pipeline_info/provenance.json` (commit, tag, command, settings, status, container sha256).
- CI on pull requests: `nextflow lint` and nf-test stub runs (`.github/workflows/ci.yml`).
- Stubs for the CellBender branch (`get_batches`, `split_adata`, `run_cellbender`, `join_adatas`,
  `plot_barcode_rank`) and a stub test with `skip_cellbender: false`.
- Restructured the pipeline into nf-core directory conventions: `main.nf`
  entry point, `workflows/qc_pipeline.nf` (named workflow), `workflows/outlier_removal.nf`
  (standalone second entry point), `modules/local/<name>/main.nf` (one process
  per module directory), `conf/base.config` + per-profile `conf/*.config`,
  `nextflow_schema.json`.
- Added `params.celltype_csv` to `workflows/outlier_removal.nf`, replacing a
  hardcoded absolute path to `test/celltype_test.csv` in the previous
  `outlier_removal.nf`, so the pipeline is portable outside its original
  checkout location.

### `Notes`

- Ported from the flat pipeline in the parent `nextflow/` directory
  (2026-08-24). Only the modules/scripts actually used by
  `qc-pipeline.nf` and `outlier_removal.nf` were carried over.
