#!/usr/bin/env nextflow

/*
 * Split adata by batch
 */
process split_adata {

    container "${params.container_dir}/utricle-qc.sif"

    input :
    path raw_adata
    val batch
    val batch_col

    output :
    tuple val(batch), path("adata_split_${batch}.h5ad"), emit : split_path

    script :
    """
    split_adata.py --adata ${raw_adata} --batch ${batch} --batch-col ${batch_col}
    """

    stub:
    """
    touch adata_split_${batch}.h5ad
    """
}
