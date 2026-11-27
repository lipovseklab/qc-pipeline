/*
 * cNMF: topometry on one subset (bin/topometry.py --cnmf_outputs). Selects the HVGs that cNMF uses and writes the
 * raw counts (genes in >= 3 cells) that cnmf prepare reads, plus the topometry AnnData for the analysis step.
 */

process cnmf_topometry {
    label 'cnmf_topometry'
    container "${params.container_dir}/utricle-qc.sif"
    tag "${meta.tag}"
    // files are named by the label (as before B3b: topometry.py also uses it to name earlier embeddings) and
    // published under the tag
    publishDir { "${params.outdir}/cnmf/${meta.tag}/anndatas" }, pattern: 'adata_topometry_*.h5ad', mode: 'copy',
        saveAs: { _fn -> "adata_topometry_${meta.tag}.h5ad" }
    publishDir { "${params.outdir}/cnmf/${meta.tag}" }, pattern: 'hvgs_*.csv', mode: 'copy',
        saveAs: { _fn -> "hvgs_${meta.tag}.csv" }
    publishDir { "${params.outdir}/cnmf/${meta.tag}/topometry" }, pattern: '{*.png,*.pkl}', mode: 'copy'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("counts_${meta.label}.h5ad"), path("hvgs_${meta.label}.csv"), emit: counts
    tuple val(meta), path("adata_topometry_${meta.label}.h5ad"),                     emit: adata
    path '*.png', optional: true
    path '*.pkl', optional: true

    script:
    """
    topometry.py \\
        --adata ${adata} \\
        --descriptor ${meta.label} \\
        ${meta.column ? "--metadata ${meta.column}" : ''} \\
        --cnmf_outputs
    """

    stub:
    """
    touch counts_${meta.label}.h5ad hvgs_${meta.label}.csv adata_topometry_${meta.label}.h5ad
    touch pca_variance_ratio_${meta.label}.png topometry_object_${meta.label}.pkl
    """
}
