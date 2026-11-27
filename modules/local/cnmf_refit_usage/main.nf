/*
 * Repeats cNMF's final usage refit with the genes in matching order and checks that cNMF's own usages agree.
 * They disagree when the run folder was prepared from gene names in pandas' `string` type: cNMF 1.7.1 then pairs
 * the spectra and the data on different genes (docs/cnmf.md). The pipeline now stores names as plain strings, so
 * for new folders both agree; the step fails otherwise. Writes
 * <run_name>.usages.k_<k>.dt_<threshold>.consensus.refit.txt next to the cNMF files, which are not changed.
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
