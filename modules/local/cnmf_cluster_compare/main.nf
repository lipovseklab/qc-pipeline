/*
 * Cluster comparison plots (one per clustering series: pca, topo, topo_ms) and the table of clusterings that
 * cnmf_clustree plots. Needs only the topometry AnnData.
 */

process cnmf_cluster_compare {
    label 'cnmf_small'
    container "${params.container_dir}/cnmf-analysis.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/cnmf/${meta.tag}/clustering" }, pattern: 'cluster_comparison_on_embedding_*.png', mode: 'copy'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path('cluster_table.csv'), path('clustree_series.txt'), emit: table
    path 'cluster_comparison_on_embedding_*.png'

    script:
    """
    # matplotlib / numba want writable cache folders; the container's home may not be
    export MPLCONFIGDIR=\$PWD/.matplotlib NUMBA_CACHE_DIR=\$PWD/.numba

    cluster_compare.py \\
        --adata ${adata} \\
        --color ${meta.cluster_color} \\
        --projection ${meta.projection}
    """

    stub:
    """
    touch cluster_table.csv clustree_series.txt cluster_comparison_on_embedding_topo_ms.png
    """
}
