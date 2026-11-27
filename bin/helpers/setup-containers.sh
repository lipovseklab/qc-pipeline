#!/usr/bin/env bash
# Build or pull the pipeline's container images (run from anywhere in a clone of this repo).
#   setup-containers.sh [--force] [--dest DIR] [image ...]
# --dest: where the .sif images go (default: the clone's containers/). Use the folder the profile's container_dir
#   points to, e.g. on Myriad: setup-containers.sh --dest ~/Scratch/containers
#   images: cellbender utricle-qc vaeda r-clustering cnmf_env cnmf-analysis clustree (default: all)
# Images that already exist are kept, so a rerun after a failure carries on where it stopped; --force rebuilds them.
# Extra options for `apptainer build` (not pull) go in APPTAINER_BUILD_OPTS, e.g. on WSL, where builds can stop with
# "failed to mount binfmt_misc filesystem": APPTAINER_BUILD_OPTS=--ignore-subuid setup-containers.sh
# The sha256 of the images each run used is recorded in its provenance.json.
set -euo pipefail
start=$PWD                          # a relative --dest is taken from where the script was started
cd "$(dirname "$0")/../.."

declare -A how=(
    [cellbender]="pull docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2"   # cell calling
    [utricle-qc]="build containers/utricle-qc.def"                              # QC, outliers, TopOMetry
    [vaeda]="build containers/vaeda.def"                                        # vaeda (doublets)
    [r-clustering]="build containers/r-clustering.def"                          # scDblFinder, scSHC
    [cnmf_env]="pull docker://quay.io/biocontainers/cnmf:1.7.1--pyhdfd78af_0"     # cNMF subset/prepare/factorize/consensus
    [cnmf-analysis]="build containers/cnmf-analysis.def"                        # cNMF GEP analysis, cluster comparison
    [clustree]="build containers/clustree.def"                                  # cNMF clustree plots
)
all=(cellbender utricle-qc vaeda r-clustering cnmf_env cnmf-analysis clustree)

force=false
dest=containers
images=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force) force=true ;;
        --dest) [[ $# -ge 2 ]] || { echo "--dest needs a folder" >&2; exit 1; }
                dest=$2; [[ $dest == /* ]] || dest=$start/$dest; shift ;;
        *) [[ -n "${how[$1]+x}" ]] || { echo "Unknown image '$1' (choose from: ${all[*]})" >&2; exit 1; }
           images+=("$1") ;;
    esac
    shift
done
[[ ${#images[@]} -gt 0 ]] || images=("${all[@]}")
mkdir -p "$dest"

for name in "${images[@]}"; do
    image=$dest/$name.sif
    if [[ -e "$image" && $force == false ]]; then
        echo "keeping $image (exists; --force to rebuild)"
        continue
    fi
    read -r verb source <<< "${how[$name]}"
    opts=()
    # APPTAINER_BUILD_OPTS may hold several options; `apptainer pull` rejects build-only ones such as --ignore-subuid
    [[ $verb == build ]] && read -r -a opts <<< "${APPTAINER_BUILD_OPTS:-}"
    apptainer "$verb" --force "${opts[@]}" "$image" "$source"
done
