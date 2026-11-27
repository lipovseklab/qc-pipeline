#!/usr/bin/env nextflow

/*
 * clustree plot of the clusterings of TopOMetry's chosen latent space, from the cluster table that combine_clusters
 * writes (the same R script draws the cNMF steps' clustree plots)
 */

process clustree {

    label 'clustree'
    container "${params.container_dir}/clustree.sif"

    input:
    tuple path(table), path(series)

    output:
    path "clustree.png", emit : plot

    script:
    """
    clustree_series.R ${table} ${series}
    mv clustree_*.png clustree.png
    """

    stub:
    """
    touch clustree.png
    """
}
