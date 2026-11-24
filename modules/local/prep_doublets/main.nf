#!/usr/bin/env nextflow

/*
 * vaeda
 */

 process prep_doublets {

    container "${params.container_dir}/utricle-qc.sif"

    input:
    path adata

    output:
    path "adata_scDblFinder_in.h5ad", emit : scDblFinder
    path "adata_vaeda_in.h5ad", emit : vaeda

    script:
    """
    prep_doublets.py --adata ${adata}
    """

    stub:
    """
    touch adata_scDblFinder_in.h5ad
    touch adata_vaeda_in.h5ad
    """
 }
