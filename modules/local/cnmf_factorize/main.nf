/*
 * cNMF factorize: one task per chunk of worker indices (1..n_jobs, `stride` at a time). cNMF takes the index
 * modulo the total workers, so 1..N covers every run. Completed iterations are skipped on a retry.
 */

process cnmf_factorize {
    label 'cnmf_factorize'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"

    input:
    tuple val(meta), val(workers)

    output:
    val meta, emit: meta

    script:
    def outdir = file("${params.outdir}/cnmf/${meta.tag}").toAbsolutePath()
    """
    for i in ${workers.join(' ')}; do
        cnmf factorize \\
            --output-dir ${outdir} \\
            --name ${meta.run_name} \\
            --worker-index \$i \\
            --total-workers ${meta.n_jobs} \\
            --skip-completed-runs
    done
    """

    stub:
    """
    echo "factorize ${meta.tag}: workers ${workers.join(' ')}"
    """
}
