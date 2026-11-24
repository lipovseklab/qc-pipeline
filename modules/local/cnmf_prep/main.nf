/*
 * cNMF prepare. cNMF writes its run folder <outdir>/cnmf/<tag>/<run_name> itself, outside the work directory:
 * factorize, combine and consensus all add to that folder.
 */

process cnmf_prep {
    label 'cnmf_small'
    container "${params.container_dir}/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/cnmf/${meta.tag}/pipeline_info" }, pattern: 'run_name.txt', mode: 'copy'

    input:
    tuple val(meta), path(counts), path(hvgs)

    output:
    val meta, emit: meta
    path 'run_name.txt'

    script:
    def outdir = file("${params.outdir}/cnmf/${meta.tag}").toAbsolutePath()   // cnmf adds /${meta.run_name}
    """
    echo ${meta.run_name} > run_name.txt

    # cnmf factorize --skip-completed-runs keeps every spectra file already in the run folder. Spectra made with
    # another seed would then be mixed into the consensus without any trace, so the seed is recorded per run folder.
    seed_file=${outdir}/${meta.run_name}/cnmf_seed.txt
    if [ -f "\$seed_file" ] && [ "\$(cat "\$seed_file")" != "${meta.seed}" ]; then
        echo "ERROR: ${outdir}/${meta.run_name} was made with seed \$(cat "\$seed_file"), not ${meta.seed}." >&2
        echo "Use another run_name or outdir for a new seed, or set seed: \$(cat "\$seed_file")." >&2
        exit 1
    fi
    if [ ! -f "\$seed_file" ] && compgen -G "${outdir}/${meta.run_name}/cnmf_tmp/*.spectra.k_*.iter_*.df.npz" > /dev/null; then
        echo "ERROR: ${outdir}/${meta.run_name} already holds factorize results but no record of their seed." >&2
        echo "Use another run_name or outdir (or write the seed used into \$seed_file if you know it)." >&2
        exit 1
    fi

    cnmf prepare \\
        --output-dir ${outdir} \\
        --name ${meta.run_name} \\
        -c ${counts} \\
        -k ${meta.k.join(' ')} \\
        --n-iter ${meta.n_iters} \\
        --genes-file ${hvgs} \\
        --seed ${meta.seed} \\
        --total-workers ${meta.n_jobs}
    echo ${meta.seed} > "\$seed_file"
    """

    stub:
    """
    echo ${meta.run_name} > run_name.txt
    """
}
