#!/usr/bin/env nextflow

/*
 * Join cellbender results to raw data
 */
process join_adatas {

    container "${params.container_dir}/utricle-qc.sif"

    input :
    path raw_adata
    path cellbender_paths

    output :
    path 'adata_cellbender_filtered.h5ad',   emit : filtered
    path 'adata_cellbender_unfiltered.h5ad', emit : all

    script :
    """
    join_adatas.py --raw-adata ${raw_adata} ${cellbender_paths}
    """
}
