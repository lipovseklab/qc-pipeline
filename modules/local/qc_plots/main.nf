#!/usr/bin/env nextflow

process qc_plots {

    container "${params.container_dir}/utricle-qc.sif"

    input:
    path adata
    val batch
    val descriptor

    output:
    path "qcplot_*${descriptor}.png"

    script:
    """
    qc_plots.py --adata ${adata} --batch ${batch} --descriptor ${descriptor}
    """

    stub:
    """
    touch qcplot_stub_${descriptor}.png
    """
}
