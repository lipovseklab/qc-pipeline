/*
* Modules
*/
include { plt_proportion }    from '../modules/local/plt_proportion/main.nf'
include { outlier_filter }    from '../modules/local/outlier_filter/main.nf'
include { filtering_summary } from '../modules/local/filtering_summary/main.nf'
include { report_stats }  from '../modules/local/report_stats/main.nf'
include { qc_plots }      from '../modules/local/qc_plots/main.nf'
include { run_topometry } from '../modules/local/run_topometry/main.nf'
include { cell_typing }   from '../modules/local/cell_typing/main.nf'
include { clustree }         from '../modules/local/clustree/main.nf'
include { scSHC }            from '../modules/local/scshc/main.nf'
include { combine_clusters } from '../modules/local/combine_clusters/main.nf'

/*
* Workflow
*/
workflow OUTLIERS {

    take:
    adata              // path: adata after QC + doublet calling (adata_doublets.h5ad)
    ds                 // dataset settings (see main.nf)
    ch_prior_stats     // channel: stats CSVs from the QC stages (report_stats_<stage>.csv)
    ch_prior_batch     // channel: per-batch stats CSVs from the QC stages (report_stats_<stage>_batch.csv)

    main:

    //////////////////////////////////////////////////
    // Outlier removal & reporting
    //////////////////////////////////////////////////

    // outlier filter (returns data in adata.raw if supplied)
    outlier_filter(
        adata,
        ds.clusters,
        ds.doublet_column,
        ds.name
    )

    // qc plots
    qc_plots(
        outlier_filter.out.adata,
        ds.batch,
        'outlierFilter'
    )

    // report stats
    report_stats(
        outlier_filter.out.adata,
        'outlierFilter',
        ds.batch
    )

    // QC-stage stats + this stage's stats, for the filtering summary across all stages
    ch_csvs = ch_prior_stats
                .mix(report_stats.out.stats)
                .collect()
    ch_batch_csvs = ch_prior_batch
                        .mix(report_stats.out.batch_stats)
                        .collect()
                        .ifEmpty([])

    // filtering_summary
    filtering_summary(
        ch_csvs,
        ch_batch_csvs
    )

    // plot final proportions per batch and cell type
    plt_proportion(
        outlier_filter.out.adata,
        ds.batch,
        ds.name
    )

    //////////////////////////////////////////////////
    // Generate the final dataset
    //////////////////////////////////////////////////

    // re-run topometry on filtered data (metadata for plots = batch column; was hard-coded 'plate')
    run_topometry(
        outlier_filter.out.adata,
        ds.batch,
        'postqc'
    )

    // scSHC clusters (optional — skip with --skip_scSHC)
    if (!ds.skip_scSHC) {
        scSHC(outlier_filter.out.adata)
    }
    scSHC_csv = ds.skip_scSHC ? channel.value('NO_FILE') : scSHC.out.csv

    // simple cell typing
    cell_typing(
        outlier_filter.out.adata,
        ds.celltype_csv,
        'postqc'
    )

    // cluster comparison plot
    combine_clusters(
        run_topometry.out.adata,
        scSHC_csv,
        cell_typing.out.csv
    )

    // clustree of the same clusterings (combine_clusters writes its input table)
    clustree(combine_clusters.out.clustree_input)

    emit:
    // outliers
    madplots                = outlier_filter.out.madplots
    csv_outliers            = outlier_filter.out.csv
    adata_outliers          = outlier_filter.out.adata
    log_outliers            = outlier_filter.out.log
    // qc reporting
    plot_qc                 = qc_plots.out
    plot_stats              = report_stats.out.stats
    plot_stats_batch        = report_stats.out.batch_stats
    // summary plots
    plot_filter_summary     = filtering_summary.out
    plot_proportions        = plt_proportion.out
    // topometry
    plot_topo_pca           = run_topometry.out.plt_pca
    plot_topo_geom          = run_topometry.out.plt_geometry
    plot_topo_embeddings    = run_topometry.out.plt_comparison
    topo_object             = run_topometry.out.tg
    // clustree
    plot_clustree           = clustree.out.plot
    // scSHC
    csv_scSHC  = ds.skip_scSHC ? channel.empty() : scSHC.out.csv
    plot_scSHC = ds.skip_scSHC ? channel.empty() : scSHC.out.plot
    rds_scSHC  = ds.skip_scSHC ? channel.empty() : scSHC.out.results_object
    // combine clusters
    plot_cluster_comparison = combine_clusters.out.plot
    adata_cluster           = combine_clusters.out.adata
}
