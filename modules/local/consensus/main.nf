// consensus writes a cache in the run folder, so combine + consensus run under an exclusive lock
process CONSENSUS {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${sample}"
    publishDir { "${params.outdir}/${sample}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(sample), val(run_name)

    output:
    tuple val(sample), val(run_name), emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file("${params.outdir}/${sample}").toAbsolutePath()   // cnmf adds /${run_name}
    """
    (
        flock -x 9
        cnmf combine --output-dir ${outdir} --name ${run_name}
        cnmf consensus \\
            --output-dir ${outdir} \\
            --name ${run_name} \\
            --components ${params.selected_k} \\
            --local-density-threshold ${params.threshold} \\
            --show-clustering
    ) 9>"${outdir}/${run_name}/.consensus.lock"
    """
}
