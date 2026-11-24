/*
 * cNMF consensus for the chosen k. consensus writes a cache in the run folder, so combine + consensus run under
 * an exclusive lock (see the local_density_cache race in the CHANGELOG).
 */

process cnmf_consensus {
    label 'cnmf_consensus'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"

    input:
    val meta

    output:
    val meta, emit: meta

    script:
    def outdir = file("${params.outdir}/cnmf/${meta.tag}").toAbsolutePath()   // cnmf adds /${meta.run_name}
    """
    (
        flock -x 9
        cnmf combine --output-dir ${outdir} --name ${meta.run_name}
        cnmf consensus \\
            --output-dir ${outdir} \\
            --name ${meta.run_name} \\
            --components ${meta.selected_k} \\
            --local-density-threshold ${meta.threshold} \\
            --show-clustering
    ) 9>"${outdir}/${meta.run_name}/.consensus.lock"
    """

    stub:
    """
    echo "consensus ${meta.tag}: k=${meta.selected_k} threshold=${meta.threshold}"
    """
}
