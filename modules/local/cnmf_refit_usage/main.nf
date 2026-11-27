/*
 * cNMF 1.7.1 writes wrong usages (its final refit pairs the spectra and the data on different genes), so the final
 * step is repeated with the genes in matching order. Reads the consensus outputs and writes
 * <run_name>.usages.k_<k>.dt_<threshold>.consensus.refit.txt next to them; the cNMF files are not changed.
 */

process cnmf_refit_usage {
    label 'cnmf_consensus'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"

    input:
    val meta

    output:
    val meta, emit: meta

    script:
    def outdir = file("${params.outdir}/cnmf/${meta.tag}").toAbsolutePath()   // the run folder is ${outdir}/${run_name}
    """
    (
        flock -x 9
        refit_usage.py \\
            --output_dir ${outdir} \\
            --name ${meta.run_name} \\
            --k ${meta.selected_k} \\
            --threshold ${meta.threshold}
    ) 9>"${outdir}/${meta.run_name}/.consensus.lock"
    """

    stub:
    """
    echo "refit usage ${meta.tag}"
    """
}
