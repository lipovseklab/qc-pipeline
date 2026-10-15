process SUBSET {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/anndatas" }, pattern: 'adata_*.h5ad', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("adata_${meta.tag}.h5ad"), emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    // the subset value may hold spaces: single-quote it for the shell; the file name uses the tag (safe label)
    def subset = "'" + meta.subset.toString().replace("'", "'\"'\"'") + "'"
    """
    subset.py \\
        --counts ${counts} \\
        --metadata ${meta.metadata} \\
        --subset ${subset} \\
        --label ${meta.label} \\
        --output-dir .

    # fail here rather than passing an empty AnnData to topometry
    python - ${meta.tag} ${meta.metadata} ${subset} ${counts} <<'PY'
    import anndata as ad, sys
    tag, metadata, subset, counts = sys.argv[1:5]
    a = ad.read_h5ad(f'adata_{tag}.h5ad', backed='r')
    print(f'subset {tag}:', a.shape)
    if a.n_obs == 0:
        sys.exit(f"Subset {metadata} == {subset} has 0 cells in {counts}")
    PY
    """
}
