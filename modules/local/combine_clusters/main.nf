#!/usr/bin/env nextflow

process combine_clusters {

    container "${params.container_dir}/utricle-qc.sif"

    input:
    path adata
    val scSHC_csv
    path celltype_csv

    output:
    path 'cluster_comparison_on_embedding.png', emit : plot
    path "adata_clusters.h5ad", emit : adata
    tuple path('cluster_table.csv'), path('clustree_series.txt'), emit : clustree_input

    script:
    """
    combine_clusters.py --adata ${adata} --celltype-csv ${celltype_csv} --scshc-csv ${scSHC_csv}
    """

    stub:
    """
    printf 'cell,topo_clusters_ms_res0.2\nc1,0\n' > cluster_table.csv
    printf 'topo_ms\ttopo_clusters_ms_res\tTopometry - multiscale latent space\n' > clustree_series.txt
    touch cluster_comparison_on_embedding.png
    touch adata_clusters.h5ad
    """
 }
