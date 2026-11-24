# utricleqc: Changelog

## [Unreleased]

### `Fixed`

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
