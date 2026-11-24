#!/usr/bin/env bash
# Build or pull the pipeline's container images into containers/ (run from anywhere in a clone of this repo).
#   setup-containers.sh [--force] [image ...]     images: cellbender utricle-qc vaeda r-clustering (default: all)
# Images that already exist are kept, so a rerun after a failure carries on where it stopped; --force rebuilds them.
# The sha256 of the images each run used is recorded in its provenance.json.
set -euo pipefail
cd "$(dirname "$0")/../.."

declare -A how=(
    [cellbender]="pull docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2"   # cell calling
    [utricle-qc]="build containers/utricle-qc.def"                              # QC, outliers, TopOMetry
    [vaeda]="build containers/vaeda.def"                                        # vaeda (doublets)
    [r-clustering]="build containers/r-clustering.def"                          # scDblFinder, scSHC, clustree
)
all=(cellbender utricle-qc vaeda r-clustering)

force=false
images=()
for arg in "$@"; do
    case "$arg" in
        --force) force=true ;;
        *) [[ -n "${how[$arg]+x}" ]] || { echo "Unknown image '$arg' (choose from: ${all[*]})" >&2; exit 1; }
           images+=("$arg") ;;
    esac
done
[[ ${#images[@]} -gt 0 ]] || images=("${all[@]}")

for name in "${images[@]}"; do
    image=containers/$name.sif
    if [[ -e "$image" && $force == false ]]; then
        echo "keeping $image (exists; --force to rebuild)"
        continue
    fi
    read -r verb source <<< "${how[$name]}"
    apptainer "$verb" --force "$image" "$source"
done
