#!/usr/bin/env bash
# Build or pull the pipeline's container images into containers/ (run from anywhere in a clone of this repo).
# Images that already exist are kept, so a rerun after a failure carries on where it stopped;
# use --force to rebuild them all. The sha256 of the images each run used is recorded in its provenance.json.
set -euo pipefail
cd "$(dirname "$0")/../.."
force=false
[[ "${1:-}" == "--force" ]] && force=true

make_image() {    # make_image <image> <pull|build> <source>
    local image=$1 how=$2 source=$3
    if [[ -e "$image" && $force == false ]]; then
        echo "keeping $image (exists; --force to rebuild)"
        return
    fi
    apptainer "$how" --force "$image" "$source"
}

make_image containers/cellbender.sif   pull  docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2   # cell calling
make_image containers/utricle-qc.sif   build containers/utricle-qc.def                                # QC, outliers, TopOMetry
make_image containers/vaeda.sif        build containers/vaeda.def                                     # vaeda (doublets)
make_image containers/r-clustering.sif build containers/r-clustering.def                              # scDblFinder, scSHC, clustree
