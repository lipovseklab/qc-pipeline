#!/usr/bin/env nextflow
/*
* lipovseklab/qc-pipeline — mouse utricle single-nucleus QC.
*
* One entry point, driven by one params file:
*   nextflow run lipovseklab/qc-pipeline -r <tag> -profile <profile> -params-file run.yaml
*
* run.yaml:
*   step: qc | outliers | all        (all = qc, then outlier removal on its output)
*       | cnmf | consensus | analysis   (cNMF on the annotated adata: three separate runs on the same outdir)
*   outdir: <folder>
*   datasets:                        (exactly one entry for now)
*     - name: alldata
*       adata: <path>                (qc: raw adata; outliers: adata_doublets.h5ad from a qc run;
*                                     cnmf: annotated adata with raw counts)
*       cnmf:                        (cnmf steps only: settings + one entry per subset, see cnmfRuns below)
*       ...                          (see nextflow_schema.json and examples/)
*/

include { validateParameters } from 'plugin/nf-schema'
include { QC }       from './workflows/qc.nf'
include { OUTLIERS } from './workflows/outliers.nf'
include { CNMF; CNMF_CONSENSUS; CNMF_ANALYSIS } from './workflows/cnmf.nf'
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

def cnmfStep() { params.step in ['cnmf', 'consensus', 'analysis'] }

def datasetSettings() {
    if( !params.datasets || params.datasets.size() != 1 )
        error "run.yaml must list exactly one dataset under `datasets:` (found ${params.datasets?.size() ?: 0})"
    def unknown = params.datasets[0].keySet() - datasetDefaults().keySet() - ['name', 'adata', 'cnmf']
    if( unknown )
        error "Unknown dataset setting(s) in run.yaml: ${unknown.join(', ')}. Allowed: name, adata, cnmf, ${datasetDefaults().keySet().join(', ')}"
    def ds = datasetDefaults() + params.datasets[0]
    // relative paths are taken from the launch folder; tasks need absolute ones
    for( k in ['adata', 'celltype_csv', 'cellbender_csv', 'qc_stats'] ) {
        if( ds[k] instanceof CharSequence && !ds[k].toString().trim().isEmpty() && !workflow.stubRun )
            ds[k] = file(ds[k].toString()).toAbsolutePath().toString()
    }
    if( ds.cnmf != null && !(ds.cnmf instanceof Map) )
        error "Dataset setting cnmf must be a block of settings (got '${ds.cnmf}')"
    // types (nf-schema can't check list entries here, see nextflow.config `validation`); the cnmf block is checked in cnmfRuns
    def wrong = ds.findAll { k, v -> k != 'cnmf' && v != null && (k.startsWith('skip_') ? !(v instanceof Boolean) : !(v instanceof CharSequence)) }
    if( wrong )
        error "Wrong type in run.yaml: " + wrong.collect { k, v -> k.startsWith('skip_') ? "${k} must be true/false (got '${v}')" : "${k} must be text (got ${v})" }.join('; ')
    // input files must exist (stub runs use placeholders)
    if( !workflow.stubRun ) {
        def absent = ['adata', 'celltype_csv', 'cellbender_csv', 'qc_stats'].findAll { ds[it] && !file(ds[it]).exists() }
        if( absent )
            error "Not found: " + absent.collect { "${it} = ${ds[it]}" }.join('; ')
    }
    def need = ['name', 'adata'] + (cnmfStep() ? [] : ['celltype_csv']) + (params.step in ['outliers', 'all'] ? ['clusters'] : [])
    def missing = need.findAll { k -> ds[k] == null || ds[k].toString().trim().isEmpty() }
    if( missing )
        error "Dataset '${ds.name ?: '?'}' is missing: ${missing.join(', ')} (step: ${params.step})"
    // the schema's patterns for dataset entries are not applied (see nextflow.config `validation`), so check them
    // here: text settings go unquoted into task commands and the name into file names
    if( !(ds.name ==~ /[A-Za-z0-9_.-]+/) )
        error "Dataset name '${ds.name}' may contain only letters, digits, '_', '.' and '-'"
    def spaced = ds.findAll { k, v -> v instanceof CharSequence && (v.toString().trim().isEmpty() || (v.toString() =~ /\s/).find()) }
    if( spaced )
        error "Dataset settings must not be blank or contain spaces: ${spaced.keySet().join(', ')}"
    if( !ds.adata.toString().endsWith('.h5ad') )
        error "adata must be an .h5ad file (got ${ds.adata})"
    if( params.step == 'outliers' && !ds.qc_stats )
        log.warn "No qc_stats folder given: the filtering summary will only show the outlier stage"
    return ds
}

/*
* cNMF settings. In run.yaml, under the dataset:
*
*   cnmf:
*     column: ctype            obs column the subsets are taken from
*     k: [20, 30, 40]          settings here apply to every subset ...
*     subsets:
*       - value: HC            ... and can be overridden per subset
*         k: [10, 15, 20]
*         selected_k: 15       consensus/analysis: the k chosen from the k selection plot
*         threshold: 0.1       consensus/analysis: local density threshold
*         clusters: topo_clusters_ms_res0.6   analysis: obs column with the clusters (from clustree)
*       - {}                   no value: the whole dataset (topometry selects the HVGs, or give `hvgs:`)
*
* Without `subsets`, the whole dataset is one run. Each run gets a folder <outdir>/cnmf/<tag>, tag = <column>_<value>
* (other characters than letters, digits, . _ - in the value become _) or the dataset name for the whole dataset.
*/
def cnmfDefaults() { [
    column        : null,           // obs column the subsets are taken from
    k             : [20, 30, 40],   // numbers of components to factorize
    n_iters       : 100,            // iterations (replicates) per k
    stride        : 10,             // factorize iterations per task
    seed          : 42,             // cnmf prepare seed
    cluster_color : 'ctype',        // obs column shown next to the clusterings (cluster comparison plot)
    selected_k    : null,           // consensus/analysis: chosen k
    threshold     : 2,              // consensus/analysis: local density threshold (2 = no filtering)
    clusters      : 'clusters',     // analysis: obs column with the cluster labels
    projection    : 'projection',   // obsm key of the embedding to plot on
    age           : 'age_pretty',   // analysis: obs column with the age labels
    summary_cols  : null,           // analysis: projection summary columns (default: age, clusters, ctype_detailed)
    usage_cutoff  : 0.1,            // analysis: usage above which a cell counts as using a GEP
    n_top_genes   : 100,            // analysis: top genes per GEP
    geps_per_page : 6,              // analysis: GEPs per page in the multi-GEP figures (0 = one page)
    formats       : 'png,pdf',      // analysis: figure formats
    analysis_dir  : 'analysis',     // analysis: output folder name in <outdir>/cnmf/<tag>
] }

def cnmfRuns(ds) {
    def c = ds.cnmf
    if( !(c instanceof Map) )
        error "step ${params.step} needs a `cnmf:` block under the dataset in run.yaml (see examples/cnmf.yaml)"
    def blockKeys  = cnmfDefaults().keySet() + ['subsets']
    def subsetKeys = cnmfDefaults().keySet() + ['value', 'run_name', 'hvgs']
    def unknown = c.keySet() - blockKeys
    if( unknown )
        error "Unknown cnmf setting(s) in run.yaml: ${unknown.join(', ')}. Allowed: ${blockKeys.join(', ')}"
    if( c.subsets != null && !(c.subsets instanceof List && c.subsets && c.subsets.every { e -> e instanceof Map }) )
        error "cnmf.subsets must be a list of entries, e.g. `- value: HC` (or leave it out for the whole dataset)"
    def entries = c.subsets ?: [[:]]
    def today = new Date().format('yyyyMMdd')

    def runs = entries.withIndex().collect { entry, i ->
        def where = "cnmf.subsets[${i}]${entry.value != null ? " (value ${entry.value})" : ''}"
        def bad = entry.keySet() - subsetKeys
        if( bad )
            error "Unknown setting(s) in ${where}: ${bad.join(', ')}. Allowed: ${subsetKeys.join(', ')}"
        def s = cnmfDefaults() + c.findAll { key, _v -> key != 'subsets' } + entry

        // numbers
        def k = s.k instanceof List ? s.k : [s.k]
        if( !k || !k.every { kv -> kv instanceof Integer && kv > 1 } )
            error "${where}: k must be a list of whole numbers above 1, e.g. [20, 30, 40] (got ${s.k})"
        def minimum = [n_iters: 1, stride: 1, n_top_genes: 1, seed: 0, geps_per_page: 0] + (s.selected_k != null ? [selected_k: 2] : [:])
        minimum.each { key, low ->
            if( !(s[key] instanceof Integer) || s[key] < low )
                error "${where}: ${key} must be a whole number of at least ${low} (got ${s[key]})"
        }
        ['threshold', 'usage_cutoff'].each { key ->
            if( !(s[key] instanceof Number) )
                error "${where}: ${key} must be a number (got ${s[key]})"
        }

        // text: goes unquoted into task commands
        def text = ['column', 'cluster_color', 'clusters', 'projection', 'age', 'formats', 'analysis_dir', 'run_name', 'hvgs']
        text.each { key ->
            def v = s[key]
            if( v != null && !(v instanceof CharSequence && v.toString() ==~ /[\w.,:\/+@%=~-]+/) )
                error "${where}: ${key} must be text without spaces or shell characters (got '${v}')"
        }
        if( !(s.analysis_dir ==~ /[\w.-]+/) )
            error "${where}: analysis_dir may only use letters, digits, _, . and - (got '${s.analysis_dir}')"
        def summary = s.summary_cols == null ? null : (s.summary_cols instanceof List ? s.summary_cols : [s.summary_cols])
        if( summary != null && !summary.every { col -> col instanceof CharSequence && col.toString() ==~ /[\w.-]+/ } )
            error "${where}: summary_cols must be a list of obs column names (got ${s.summary_cols})"

        // what to run on
        def value = entry.value
        if( value != null && !(value instanceof CharSequence || value instanceof Number) )
            error "${where}: value must be a single value of the ${s.column ?: 'cnmf.column'} column"
        if( value != null ) {
            def v = value.toString()
            if( v.isEmpty() || v.trim() != v || (v =~ /[\t\r\n]/).find() )
                error "${where}: value may not be empty, start or end with spaces, or hold tabs or line breaks"
        }
        if( value != null && !s.column )
            error "${where}: a subset needs `column` (in the cnmf block or the entry)"
        if( value != null && s.hvgs )
            error "${where}: hvgs is only used for the whole dataset (subsets get their HVGs from topometry)"
        if( s.hvgs && !workflow.stubRun && !file(s.hvgs).exists() )
            error "${where}: hvgs file not found: ${s.hvgs}"
        def label = value != null ? value.toString().replaceAll(/[^\w.-]+/, '_') : ds.name
        def tag   = value != null ? "${s.column}_${label}" : ds.name

        // run name: as given, else the one recorded by an earlier cnmf run into this outdir, else a new one
        def recorded = file("${params.outdir}/cnmf/${tag}/pipeline_info/run_name.txt")
        def run_name = s.run_name ?: (recorded.exists() ? recorded.text.trim() : null)
        if( !run_name ) {
            if( params.step == 'cnmf' || workflow.stubRun ) run_name = "cnmf_${tag}_${today}"
            else if( s.selected_k == null ) run_name = null      // skipped by consensus/analysis below
            else error "${where}: no cNMF run found for '${tag}' (${recorded} is missing): run step cnmf into this outdir first, or give run_name"
        }
        if( run_name != null && !(run_name ==~ /[\w.-]+/) )
            error "${where}: run name '${run_name}' may only use letters, digits, _, . and -"

        def n_jobs = k.size() * s.n_iters
        s + [
            tag          : tag,
            value        : value,
            label        : label,
            k            : k,
            hvgs         : s.hvgs,
            run_name     : run_name,
            summary_cols : summary,
            n_jobs       : n_jobs,
            n_chunks     : Math.ceil(n_jobs / s.stride) as int,
        ]
    }

    def dupTag = runs.countBy { r -> r.tag }.findAll { _t, n -> n > 1 }.keySet()
    if( dupTag )
        error "cnmf.subsets give the same folder name more than once (values that differ only in spaces or punctuation?): ${dupTag.join(', ')}"
    def dupName = runs.findAll { r -> r.run_name != null }.countBy { r -> r.run_name }.findAll { _t, n -> n > 1 }.keySet()
    if( dupName )
        error "cnmf.subsets give the same run_name more than once: ${dupName.join(', ')}"

    if( params.step == 'cnmf' )
        return runs

    // consensus / analysis: the subsets whose k has been chosen
    def chosen = runs.findAll { r -> r.selected_k != null }
    if( !chosen )
        error "step ${params.step}: no subset has a selected_k yet. Add `selected_k:` (chosen from the k selection plot) to the subsets to run"
    runs.findAll { r -> r.selected_k == null }.each { r -> log.info "cnmf: skipping ${r.tag} (no selected_k)" }
    chosen.each { m ->
        if( !(m.selected_k in m.k) )
            log.warn "cnmf ${m.tag}: selected_k ${m.selected_k} is not one of k ${m.k}"
        def run_dir = file("${params.outdir}/cnmf/${m.tag}/${m.run_name}")
        if( !workflow.stubRun && !run_dir.exists() )
            error "cnmf ${m.tag}: run folder ${run_dir} not found (run step cnmf into this outdir first)"
    }
    if( params.step == 'analysis' ) {
        chosen = chosen.collect { m ->
            // subsets and topometry runs use the topometry AnnData; a run with a given HVG list uses the dataset's adata
            def a = m.value != null || m.hvgs == null
                ? "${params.outdir}/cnmf/${m.tag}/anndatas/adata_topometry_${m.tag}.h5ad"
                : ds.adata
            if( !workflow.stubRun && !file(a).exists() )
                error "cnmf ${m.tag}: ${a} not found (run step cnmf into this outdir first)"
            m + [analysis_adata: a]
        }
        ['gmt', 'collectri'].each { key ->
            if( !workflow.stubRun && !file(params[key]).exists() )
                error "The analysis step needs ${key} = ${params[key]}: run `python bin/fetch_resources.py` once (needs internet), or set ${key}"
        }
    }
    return chosen
}

workflow {

    main:
    // Nextflow ignores a stray word without complaint, so a typo such as `threshold 0.18` (no --) would run with the
    // default value. The same happens to `-resume <run name>`, which only accepts `last` or a full session id.
    // (Not checked in stub runs: Nextflow 25.10 passes `-stub-run` on as a stray `true`.)
    if( args && !workflow.stubRun )
        error "Unexpected argument(s): ${args.join(' ')}\n" +
              "Pipeline options need two dashes (--outdir x). Nextflow options take one, and -resume only accepts a full session id (-resume 1f50da21-...), not a run name."
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
    def runs     = cnmfStep() ? cnmfRuns(ds) : []

    if( params.step == 'cnmf' ) {
        CNMF(runs, file(ds.adata))
    }
    if( params.step == 'consensus' ) {
        CNMF_CONSENSUS(runs)
    }
    if( params.step == 'analysis' ) {
        CNMF_ANALYSIS(runs, file(params.gmt), file(params.collectri))
    }

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
