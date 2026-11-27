/*
 * cNMF combine + k selection plot (<outdir>/cnmf/<tag>/<run_name>/<run_name>.k_selection.png).
 */

process cnmf_combine {
    label 'cnmf_combine'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"

    input:
    val meta

    output:
    val meta, emit: meta

    script:
    def outdir = file("${params.outdir}/cnmf/${meta.tag}").toAbsolutePath()
    """
    cnmf combine --output-dir ${outdir} --name ${meta.run_name}
    cnmf k_selection_plot --output-dir ${outdir} --name ${meta.run_name}
    """

    stub:
    """
    echo "combine ${meta.tag}"
    """
}
