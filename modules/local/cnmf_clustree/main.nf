/*
 * Cluster stability across resolutions (clustree), one plot per clustering series, from the table written by
 * cnmf_cluster_compare.
 */

process cnmf_clustree {
    label 'cnmf_small'
    container "${params.container_dir}/clustree.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/cnmf/${meta.tag}/clustering" }, pattern: 'clustree_*.png', mode: 'copy'

    input:
    tuple val(meta), path(table), path(series)

    output:
    path 'clustree_*.png'

    script:
    """
    clustree_series.R ${table} ${series}
    """

    stub:
    """
    touch clustree_topo_ms.png
    """
}
