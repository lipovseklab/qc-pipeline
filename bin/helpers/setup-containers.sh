#!/usr/bin/env bash
# Build or pull the pipeline's container images into containers/ (run from the pipeline folder).
# The paper runs used images built on 2026-03-31; their sha256 are recorded in each run's provenance.json.
set -euo pipefail
cd "$(dirname "$0")/../.."

# CellBender 0.3.2 (cell calling)
apptainer pull containers/cellbender.sif docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2

# QC, outliers, TopOMetry
apptainer build containers/utricle-qc.sif containers/utricle-qc.def

# vaeda (doublets)
apptainer build containers/vaeda.sif containers/vaeda.def

# R: scDblFinder, scSHC, clustree
apptainer build containers/r-clustering.sif containers/r-clustering.def
