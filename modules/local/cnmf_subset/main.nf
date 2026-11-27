/*
 * cNMF: keep the cells of one subset (obs[column] == value). Fails if the subset has 0 cells.
 */

process cnmf_subset {
    label 'cnmf_small'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/cnmf/${meta.tag}/anndatas" }, pattern: 'adata_*.h5ad', mode: 'copy'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("adata_${meta.tag}.h5ad"), emit: adata

    script:
    // the subset value may hold spaces: single-quote it for the shell; the file name uses the tag (safe label)
    def value = "'" + meta.value.toString().replace("'", "'\"'\"'") + "'"
    """
    subset.py \\
        --counts ${adata} \\
        --metadata ${meta.column} \\
        --subset ${value} \\
        --label ${meta.label} \\
        --output-dir .

    # fail here rather than passing an empty AnnData to topometry
    python - ${meta.tag} ${meta.column} ${value} ${adata} <<'PY'
    import anndata as ad, sys
    tag, column, value, adata = sys.argv[1:5]
    a = ad.read_h5ad(f'adata_{tag}.h5ad', backed='r')
    print(f'subset {tag}:', a.shape)
    if a.n_obs == 0:
        sys.exit(f"Subset {column} == {value} has 0 cells in {adata}")
    PY
    """

    stub:
    """
    touch adata_${meta.tag}.h5ad
    """
}
