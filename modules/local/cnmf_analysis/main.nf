/*
 * GEP analysis of one subset after consensus: merges the GEP usages into the AnnData, top genes, GO BP + CollecTRI
 * enrichment, DGE by cluster and the figures. Writes <outdir>/cnmf/<tag>/<analysis_dir>/.
 */

process cnmf_analysis {
    label 'cnmf_analysis'
    container "${params.container_dir}/cnmf-analysis.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/cnmf/${meta.tag}" }, mode: 'copy'

    input:
    tuple val(meta), path(adata)
    path gmt
    path collectri

    output:
    path "${meta.analysis_dir}"

    script:
    def run_dir = file("${params.outdir}/cnmf/${meta.tag}/${meta.run_name}").toAbsolutePath()
    """
    # matplotlib / numba want writable cache folders; the container's home may not be
    export MPLCONFIGDIR=\$PWD/.matplotlib NUMBA_CACHE_DIR=\$PWD/.numba

    # what started this run, for the call record written next to the outputs
    cat > nextflow_info.txt <<'NXF_INFO'
command: ${workflow.commandLine}
run name: ${workflow.runName}
session id: ${workflow.sessionId}
script revision: ${workflow.scriptId}
launch directory: ${workflow.launchDir}
project directory: ${workflow.projectDir}
NXF_INFO

    gep_analysis.py \\
        --adata ${adata} \\
        --run_dir ${run_dir} \\
        --run_name ${meta.run_name} \\
        --k ${meta.selected_k} \\
        --threshold ${meta.threshold} \\
        --gmt ${gmt} \\
        --collectri ${collectri} \\
        --sample ${meta.tag} \\
        --out ${meta.analysis_dir} \\
        --projection ${meta.projection} \\
        --clusters ${meta.clusters} \\
        --age ${meta.age} \\
        ${meta.summary_cols ? "--summary_cols ${meta.summary_cols.join(' ')}" : ''} \\
        --usage_cutoff ${meta.usage_cutoff} \\
        --n_top_genes ${meta.n_top_genes} \\
        --geps_per_page ${meta.geps_per_page} \\
        --formats ${meta.formats} \\
        --nextflow_info nextflow_info.txt
    """

    stub:
    """
    mkdir -p ${meta.analysis_dir}/figures ${meta.analysis_dir}/csvs
    touch ${meta.analysis_dir}/csvs/top${meta.n_top_genes}_genes.csv ${meta.analysis_dir}/analysis_call.txt
    """
}
