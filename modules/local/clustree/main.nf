#!/usr/bin/env nextflow

/*
 * clustree plot
 */

 process clustree {

    container "${params.container_dir}/r-clustering.sif"

    input:
    path adata
    path celltype_csv

    output:
    path "clustree.png", emit : plot

    script:
    """
    clustree.R ${adata} ${celltype_csv}
    """

    stub:
    """
    touch clustree.png
    """
 }
