process PREP {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/pipeline_info" }, pattern: 'run_name.txt', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(counts), path(hvgs)

    output:
    val meta, emit: res
    path 'run_name.txt'
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file("${params.outdir}/${meta.tag}").toAbsolutePath()   // cnmf adds /${meta.run_name}
    """
    echo ${meta.run_name} > run_name.txt

    # cnmf factorize --skip-completed-runs keeps every spectra file already in the run folder. Spectra made with
    # another seed would then be mixed into the consensus without any trace, so the seed is recorded per run folder.
    seed_file=${outdir}/${meta.run_name}/cnmf_seed.txt
    if [ -f "\$seed_file" ] && [ "\$(cat "\$seed_file")" != "${params.seed}" ]; then
        echo "ERROR: ${outdir}/${meta.run_name} was made with seed \$(cat "\$seed_file"), not ${params.seed}." >&2
        echo "Use another run_name or --outdir for a new seed, or rerun with --seed \$(cat "\$seed_file")." >&2
        exit 1
    fi
    if [ ! -f "\$seed_file" ] && compgen -G "${outdir}/${meta.run_name}/cnmf_tmp/*.spectra.k_*.iter_*.df.npz" > /dev/null; then
        echo "ERROR: ${outdir}/${meta.run_name} already holds factorize results but no record of their seed." >&2
        echo "Use another run_name or --outdir (or write the seed used into \$seed_file if you know it)." >&2
        exit 1
    fi

    cnmf prepare \\
        --output-dir ${outdir} \\
        --name ${meta.run_name} \\
        -c ${counts} \\
        -k ${meta.k_vals.join(' ')} \\
        --n-iter ${meta.n_iters} \\
        --genes-file ${hvgs} \\
        --seed ${params.seed} \\
        --total-workers ${meta.n_jobs}
    echo ${params.seed} > "\$seed_file"
    """
}
