/*
* cNMF on one dataset, one chain per entry under `cnmf.subsets` in run.yaml (see main.nf `cnmfRuns`).
*
*   step cnmf       subset -> topometry -> prepare -> factorize -> combine (+ k selection plot), plus clustree
*   step consensus  consensus + usage refit, for the subsets with a `selected_k`
*   step analysis   GEP analysis, for the subsets with a `selected_k`
*
* cNMF keeps its run folders under <outdir>/cnmf/<tag>/<run_name>, outside the Nextflow work directory, so the
* three steps are separate runs on the same outdir.
*/

include { cnmf_subset          } from '../modules/local/cnmf_subset/main.nf'
include { cnmf_topometry       } from '../modules/local/cnmf_topometry/main.nf'
include { cnmf_prep            } from '../modules/local/cnmf_prep/main.nf'
include { cnmf_factorize       } from '../modules/local/cnmf_factorize/main.nf'
include { cnmf_combine         } from '../modules/local/cnmf_combine/main.nf'
include { cnmf_consensus       } from '../modules/local/cnmf_consensus/main.nf'
include { cnmf_refit_usage     } from '../modules/local/cnmf_refit_usage/main.nf'
include { cnmf_analysis        } from '../modules/local/cnmf_analysis/main.nf'
include { cnmf_cluster_compare } from '../modules/local/cnmf_cluster_compare/main.nf'
include { cnmf_clustree        } from '../modules/local/cnmf_clustree/main.nf'

// subset and topometry only get the fields they use: Nextflow hashes val inputs, so the whole run settings
// (k, n_iters, stride, ...) would otherwise invalidate their cache whenever one of those changed
def smallMeta(m) {
    [tag: m.tag, column: m.column, value: m.value, label: m.label]
}

workflow CNMF {

    take:
    runs     // list of run settings (main.nf cnmfRuns)
    adata    // path: the dataset's AnnData (raw counts)

    main:
    ch_runs = channel.fromList(runs)

    ch_input = ch_runs.branch { m ->
        subset    : m.value != null                 // subset, then topometry
        topometry : m.hvgs == null                  // whole dataset, topometry for the HVGs
        direct    : true                            // whole dataset with a given HVG list: straight to prepare
    }

    cnmf_subset( ch_input.subset.map { m -> tuple(smallMeta(m), adata) } )
    cnmf_topometry( cnmf_subset.out.adata.mix( ch_input.topometry.map { m -> tuple(smallMeta(m), adata) } ) )

    // the full settings are joined back on the tag
    ch_by_tag = ch_runs.map { m -> tuple(m.tag, m) }
    ch_prep = cnmf_topometry.out.counts
        .map { sm, counts, hvgs -> tuple(sm.tag, counts, hvgs) }
        .join(ch_by_tag)
        .map { _tag, counts, hvgs, m -> tuple(m, counts, hvgs) }
        .mix( ch_input.direct.map { m -> tuple(m, adata, file(m.hvgs)) } )
    cnmf_prep(ch_prep)

    // chunks of `stride` worker indices per task
    cnmf_factorize( cnmf_prep.out.meta.flatMap { m ->
        (1..m.n_jobs).collate(m.stride).collect { chunk -> tuple(m, chunk) }
    } )

    // each run's combine starts as soon as its own chunks are done
    cnmf_combine( cnmf_factorize.out.meta
        .map { m -> tuple(groupKey(m.tag, m.n_chunks), m) }
        .groupTuple()
        .map { _key, ms -> ms[0] } )

    // cluster comparison and clustree only need the topometry AnnData, so they start as soon as it exists
    ch_topo = cnmf_topometry.out.adata
        .map { sm, a -> tuple(sm.tag, a) }
        .join(ch_by_tag)
        .map { _tag, a, m -> tuple(m, a) }
    cnmf_cluster_compare(ch_topo)
    cnmf_clustree(cnmf_cluster_compare.out.table)

    emit:
    combined = cnmf_combine.out.meta
}

workflow CNMF_CONSENSUS {

    take:
    runs     // runs with a selected_k

    main:
    cnmf_consensus( channel.fromList(runs) )
    cnmf_refit_usage(cnmf_consensus.out.meta)    // cNMF 1.7.1's own usages are wrong, see the module

    emit:
    refit = cnmf_refit_usage.out.meta
}

workflow CNMF_ANALYSIS {

    take:
    runs       // runs with a selected_k, each with `analysis_adata` set (main.nf)
    gmt        // path: GO BP gene sets (bin/fetch_resources.py)
    collectri  // path: CollecTRI network (bin/fetch_resources.py)

    main:
    cnmf_analysis(
        channel.fromList(runs).map { m -> tuple(m, file(m.analysis_adata)) },
        gmt,
        collectri
    )

    emit:
    analysis = cnmf_analysis.out
}
