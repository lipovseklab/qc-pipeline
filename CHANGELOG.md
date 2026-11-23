# utricleqc: Changelog

## [Unreleased]

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
