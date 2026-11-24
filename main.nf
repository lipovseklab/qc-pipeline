#!/usr/bin/env nextflow
/*
* lipovseklab/qc-pipeline — mouse utricle single-nucleus QC.
*
* One entry point, driven by one params file:
*   nextflow run lipovseklab/qc-pipeline -r <tag> -profile <profile> -params-file run.yaml
*
* run.yaml:
*   step: qc | outliers | all        (all = qc, then outlier removal on its output)
*   outdir: <folder>
*   datasets:                        (exactly one entry for now)
*     - name: alldata
*       adata: <path>                (qc: raw adata; outliers: adata_doublets.h5ad from a qc run)
*       ...                          (see nextflow_schema.json and examples/)
*/

include { validateParameters } from 'plugin/nf-schema'
include { QC }       from './workflows/qc.nf'
include { OUTLIERS } from './workflows/outliers.nf'
include { checkCodeState; writeProvenance; containerHashes } from './subworkflows/local/provenance/main.nf'

// Settings every dataset gets unless run.yaml says otherwise. Keep in sync with nextflow_schema.json.
def datasetDefaults() { [
    batch           : 'plate',              // obs column with the batch (plate)
    sample          : 'sample',             // obs column used for QC plots
    celltype_csv    : null,                 // marker table for simple cell typing (required)
    skip_cellbender : true,                 // true when the adata is already cell-called
    cellbender_csv  : null,                 // optional per-batch CellBender settings
    skip_scSHC      : false,                // scSHC is memory hungry
    clusters        : null,                 // outliers: obs column with the clusters to compute MADs in
    doublet_column  : 'scDblFinder_calls',  // outliers: obs column with doublet calls ('doublet'/'singlet')
    qc_stats        : null,                 // outliers alone: folder with report_stats_*.csv from the qc run
] }

def datasetSettings() {
    if( !params.datasets || params.datasets.size() != 1 )
        error "run.yaml must list exactly one dataset under `datasets:` (found ${params.datasets?.size() ?: 0})"
    def unknown = params.datasets[0].keySet() - datasetDefaults().keySet() - ['name', 'adata']
    if( unknown )
        error "Unknown dataset setting(s) in run.yaml: ${unknown.join(', ')}. Allowed: name, adata, ${datasetDefaults().keySet().join(', ')}"
    def ds = datasetDefaults() + params.datasets[0]
    // types (nf-schema can't check list entries here, see nextflow.config `validation`)
    def wrong = ds.findAll { k, v -> v != null && (k.startsWith('skip_') ? !(v instanceof Boolean) : !(v instanceof CharSequence)) }
    if( wrong )
        error "Wrong type in run.yaml: " + wrong.collect { k, v -> k.startsWith('skip_') ? "${k} must be true/false (got '${v}')" : "${k} must be text (got ${v})" }.join('; ')
    // input files must exist (stub runs use placeholders)
    if( !workflow.stubRun ) {
        def absent = ['adata', 'celltype_csv', 'cellbender_csv', 'qc_stats'].findAll { ds[it] && !file(ds[it]).exists() }
        if( absent )
            error "Not found: " + absent.collect { "${it} = ${ds[it]}" }.join('; ')
    }
    def need = ['name', 'adata', 'celltype_csv'] + (params.step in ['outliers', 'all'] ? ['clusters'] : [])
    def missing = need.findAll { k -> ds[k] == null }
    if( missing )
        error "Dataset '${ds.name ?: '?'}' is missing: ${missing.join(', ')} (step: ${params.step})"
    if( params.step == 'outliers' && !ds.qc_stats )
        log.warn "No qc_stats folder given: the filtering summary will only show the outlier stage"
    return ds
}

workflow {

    main:
    validateParameters()
    def code = checkCodeState()     // stops on uncommitted pipeline code unless --allow_dirty
    def ds = datasetSettings()
    writeProvenance(code, ds, [status: 'running'])
    def wf = workflow                 // `workflow` and `params` are not visible inside the handler closure
    def container_dir = params.container_dir
    wf.onComplete {
        writeProvenance(code, ds, [
            status    : wf.success ? 'completed' : 'failed',
            completed : wf.complete?.toString(),
            duration  : wf.duration?.toString(),
            error     : wf.errorMessage,
            containers: [dir: container_dir, sha256: containerHashes()],
        ])
    }
    def qc       = params.step in ['qc', 'all']
    def outliers = params.step in ['outliers', 'all']

    if( qc ) {
        QC(ds)
    }
    if( outliers ) {
        def stats_dir = ds.qc_stats ? file(ds.qc_stats) : null
        OUTLIERS(
            qc ? QC.out.adata_doublets : file(ds.adata),
            ds,
            qc ? QC.out.stats_all
               : (stats_dir ? channel.fromPath("${stats_dir}/report_stats_*.csv").filter { !it.name.endsWith('_batch.csv') } : channel.empty()),
            qc ? QC.out.batch_stats_all
               : (stats_dir ? channel.fromPath("${stats_dir}/report_stats_*_batch.csv") : channel.empty())
        )
    }

    publish:
    qc_plot_barcode_rank         = qc ? QC.out.plot_barcode_rank : channel.empty()
    qc_adata_joined_batches      = qc ? QC.out.adata_joined_batches : channel.empty()
    qc_cellbender_h5             = qc ? QC.out.cellbender_h5 : channel.empty()
    qc_cellbender_report         = qc ? QC.out.cellbender_report : channel.empty()
    qc_cellbender_log            = qc ? QC.out.cellbender_log : channel.empty()
    qc_cellbender_pdf            = qc ? QC.out.cellbender_pdf : channel.empty()
    qc_cellbender_metrics        = qc ? QC.out.cellbender_metrics : channel.empty()
    qc_adata_doublets            = qc ? QC.out.adata_doublets : channel.empty()
    qc_csv_vaeda                 = qc ? QC.out.csv_vaeda : channel.empty()
    qc_csv_scdblfinder           = qc ? QC.out.csv_scdblfinder : channel.empty()
    qc_plot_doublets_comparison  = qc ? QC.out.plot_doublets_comparison : channel.empty()
    qc_sessionInfo_scdblfinder   = qc ? QC.out.sessionInfo_scdblfinder : channel.empty()
    qc_plot_topo_pca             = qc ? QC.out.plot_topo_pca : channel.empty()
    qc_plot_topo_geom            = qc ? QC.out.plot_topo_geom : channel.empty()
    qc_plot_topo_embeddings      = qc ? QC.out.plot_topo_embeddings : channel.empty()
    qc_plot_clustree             = qc ? QC.out.plot_clustree : channel.empty()
    qc_csv_scSHC                 = qc ? QC.out.csv_scSHC : channel.empty()
    qc_plot_scSHC                = qc ? QC.out.plot_scSHC : channel.empty()
    qc_rds_scSHC                 = qc ? QC.out.rds_scSHC : channel.empty()
    qc_plot_cluster_comparison   = qc ? QC.out.plot_cluster_comparison : channel.empty()
    qc_qcplots_unfiltered        = qc ? QC.out.qcplots_unfiltered : channel.empty()
    qc_qcplots_cellcalling       = qc ? QC.out.qcplots_cellcalling : channel.empty()
    qc_qcplots_minFilter         = qc ? QC.out.qcplots_minFilter : channel.empty()
    qc_stats_unfiltered          = qc ? QC.out.stats_unfiltered : channel.empty()
    qc_stats_cellCalling         = qc ? QC.out.stats_cellCalling : channel.empty()
    qc_stats_minFilter           = qc ? QC.out.stats_minFilter : channel.empty()
    qc_batch_stats_unfiltered    = qc ? QC.out.batch_stats_unfiltered : channel.empty()
    qc_batch_stats_cellCalling   = qc ? QC.out.batch_stats_cellCalling : channel.empty()
    qc_batch_stats_minFilter     = qc ? QC.out.batch_stats_minFilter : channel.empty()

    outliers_madplots            = outliers ? OUTLIERS.out.madplots : channel.empty()
    outliers_csv_outliers        = outliers ? OUTLIERS.out.csv_outliers : channel.empty()
    outliers_adata_outliers      = outliers ? OUTLIERS.out.adata_outliers : channel.empty()
    outliers_log_outliers        = outliers ? OUTLIERS.out.log_outliers : channel.empty()
    outliers_plot_qc             = outliers ? OUTLIERS.out.plot_qc : channel.empty()
    outliers_plot_stats          = outliers ? OUTLIERS.out.plot_stats : channel.empty()
    outliers_plot_stats_batch    = outliers ? OUTLIERS.out.plot_stats_batch : channel.empty()
    outliers_plot_filter_summary = outliers ? OUTLIERS.out.plot_filter_summary : channel.empty()
    outliers_plot_proportions    = outliers ? OUTLIERS.out.plot_proportions : channel.empty()
    outliers_plot_topo_pca       = outliers ? OUTLIERS.out.plot_topo_pca : channel.empty()
    outliers_plot_topo_geom      = outliers ? OUTLIERS.out.plot_topo_geom : channel.empty()
    outliers_plot_topo_embeddings= outliers ? OUTLIERS.out.plot_topo_embeddings : channel.empty()
    outliers_topo_object         = outliers ? OUTLIERS.out.topo_object : channel.empty()
    outliers_plot_clustree       = outliers ? OUTLIERS.out.plot_clustree : channel.empty()
    outliers_csv_scSHC           = outliers ? OUTLIERS.out.csv_scSHC : channel.empty()
    outliers_plot_scSHC          = outliers ? OUTLIERS.out.plot_scSHC : channel.empty()
    outliers_rds_scSHC           = outliers ? OUTLIERS.out.rds_scSHC : channel.empty()
    outliers_plot_cluster_comparison= outliers ? OUTLIERS.out.plot_cluster_comparison : channel.empty()
    outliers_adata_cluster       = outliers ? OUTLIERS.out.adata_cluster : channel.empty()
}

output {


    /*
    * NB: Final adata for outliers etc is combine_doublets.out.adata
    */

    // adatas
    // raw - cell calling labelled
    qc_adata_joined_batches {
        path 'qc/anndatas'
        mode 'copy'
    }
    // final output
    qc_adata_doublets {
        path 'qc/anndatas'
        mode 'copy'
    }

    // plots
    qc_plot_barcode_rank {
        path 'qc/plots/cellcalling'
        mode 'copy'
    }
    qc_plot_doublets_comparison {
        path 'qc/plots/doublets'
        mode 'copy'
    }
    qc_plot_cluster_comparison {
        path 'qc/plots/clusters'
        mode 'copy'
    }
    qc_plot_clustree {
        path 'qc/plots/clusters'
        mode 'copy'
    }
    qc_plot_scSHC {
        path 'qc/plots/clusters'
        mode 'copy'
    }

    // cellbender
    qc_cellbender_h5 {
        path 'qc/anndatas/cellbender'
        mode 'copy'
    }
    qc_cellbender_report {
        path 'qc/reports/cellbender'
        mode 'copy'
    }
    qc_cellbender_log {
        path 'qc/reports/cellbender'
        mode 'copy'
    }
    qc_cellbender_pdf {
        path 'qc/reports/cellbender'
        mode 'copy'
    }
    qc_cellbender_metrics {
        path 'qc/reports/cellbender'
        mode 'copy'
    }

    // filter reports
    qc_stats_unfiltered {
        path 'qc/reports/stats'
        mode 'copy'
    }
    qc_stats_cellCalling {
        path 'qc/reports/stats'
        mode 'copy'
    }
    qc_stats_minFilter {
        path 'qc/reports/stats'
        mode 'copy'
    }
    qc_batch_stats_unfiltered {
        path 'qc/reports/stats'
        mode 'copy'
    }
    qc_batch_stats_cellCalling {
        path 'qc/reports/stats'
        mode 'copy'
    }
    qc_batch_stats_minFilter {
        path 'qc/reports/stats'
        mode 'copy'
    }

    // other reports
    qc_sessionInfo_scdblfinder {
        path 'qc/reports/doublets'
        mode 'copy'
    }
    qc_rds_scSHC {
        path 'qc/reports/clusters'
        mode 'copy'
    }

    // qc plots
    qc_qcplots_unfiltered {
        path 'qc/plots/filtering'
        mode 'copy'
    }
    qc_qcplots_cellcalling {
        path 'qc/plots/filtering'
        mode 'copy'
    }
    qc_qcplots_minFilter {
        path 'qc/plots/filtering'
        mode 'copy'
    }

    // topometry
    qc_plot_topo_pca {
        path 'qc/plots/topometry'
        mode 'copy'
    }
    qc_plot_topo_geom {
        path 'qc/plots/topometry'
        mode 'copy'
    }
    qc_plot_topo_embeddings {
        path 'qc/plots/topometry'
        mode 'copy'
    }

    // csvs
    qc_csv_vaeda {
        path 'qc/csv'
        mode 'copy'
    }
    qc_csv_scdblfinder {
        path 'qc/csv'
        mode 'copy'
    }
    qc_csv_scSHC {
        path 'qc/csv'
        mode 'copy'
    }

    // ---- outlier removal ----
    outliers_madplots {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_csv_outliers {
        path 'outliers/reports'
        mode 'copy'
    }
    outliers_adata_outliers {
        path 'outliers/anndatas'
        mode 'copy'
    }
    outliers_log_outliers {
        path 'outliers/reports'
        mode 'copy'
    }
    outliers_plot_qc {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_plot_stats {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_plot_stats_batch {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_plot_filter_summary {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_plot_proportions {
        path 'outliers/plots/qc'
        mode 'copy'
    }
    outliers_plot_topo_pca {
        path 'outliers/plots/topo'
        mode 'copy'
    }
    outliers_plot_topo_geom {
        path 'outliers/plots/topo'
        mode 'copy'
    }
    outliers_plot_topo_embeddings {
        path 'outliers/plots/topo'
        mode 'copy'
    }
    outliers_topo_object {
        path 'outliers/anndatas'
        mode 'copy'
    }
    outliers_plot_clustree {
        path 'outliers/plots/clusters'
        mode 'copy'
    }
    outliers_csv_scSHC {
        path 'outliers/reports'
        mode 'copy'
    }
    outliers_plot_scSHC {
        path 'outliers/plots/clusters'
        mode 'copy'
    }
    outliers_rds_scSHC {
        path 'outliers/reports'
        mode 'copy'
    }
    outliers_plot_cluster_comparison {
        path 'outliers/plots/clusters'
        mode 'copy'
    }
    outliers_adata_cluster {
        path 'outliers/anndatas'
        mode 'copy'
    }
}
