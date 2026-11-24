/*
* Provenance: refuse to run uncommitted code, and record exactly what ran.
*
*   checkCodeState()           call first in the entry workflow; stops on uncommitted changes unless --allow_dirty
*   writeProvenance(...)       writes <outdir>/pipeline_info/provenance.json (at the start, and again on completion)
*
* With `nextflow run lipovseklab/qc-pipeline -r <tag>` the code is a clean checkout of that tag, so the check passes and
* provenance.json names the tag and commit. Running a local folder with edits needs --allow_dirty, and the edited
* files are listed in provenance.json.
*/

def git(List args) {
    try {
        def proc = (['git', '-C', projectDir.toString()] + args).execute()
        def out = new StringBuilder()
        def err = new StringBuilder()
        proc.waitForProcessOutput(out, err)
        return proc.exitValue() == 0 ? out.toString().trim() : null
    } catch( Exception _e ) {
        return null   // git not installed
    }
}

// commit, tag and uncommitted (tracked) files of the pipeline code.
// `unknown` is set when git can't tell whether the code is clean: then the results can't be tied to a commit either.
def codeState() {
    def commit = git(['rev-parse', 'HEAD'])
    def unknown = null
    // a copy of the pipeline unpacked inside some other repository would resolve to that repository's commit
    if( commit && git(['ls-files', '--error-unmatch', '--', workflow.scriptFile.name]) == null ) {
        unknown = "${workflow.scriptFile.name} is not tracked by the git repository that holds ${projectDir}"
        commit = null
    }
    def status = commit ? git(['status', '--porcelain', '--untracked-files=no']) : null
    if( commit && status == null )
        unknown = "`git status` failed in ${projectDir}"
    [
        commit   : commit ?: workflow.commitId,
        revision : workflow.revision ?: (commit ? git(['rev-parse', '--abbrev-ref', 'HEAD']) : null),
        describe : commit ? git(['describe', '--tags', '--always', '--dirty']) : null,
        remote   : workflow.repository ?: (commit ? git(['config', '--get', 'remote.origin.url']) : null),
        changed  : status ? status.readLines()*.trim().findAll() : [],
        unknown  : unknown,
    ]
}

def checkCodeState() {
    def code = codeState()
    // clean means: a commit that git itself confirmed (workflow.commitId alone is not checked against the files)
    def confirmed = code.commit && !code.unknown && code.describe != null
    if( workflow.stubRun || params.allow_dirty ) {
        if( params.allow_dirty && !workflow.stubRun && (code.changed || !confirmed) )
            log.warn "Running with --allow_dirty: these results can't be tied to a commit (${code.unknown ?: 'edited files are listed in provenance.json'})"
        return code
    }
    if( code.unknown )
        error "The state of the pipeline code can't be checked: ${code.unknown}.\n" +
              "Run from git (`nextflow run lipovseklab/qc-pipeline -r <tag>`), or add --allow_dirty true for a test run."
    if( !confirmed )
        error "The pipeline code in ${projectDir} is not a git checkout, so these results could not be tied to a commit.\n" +
              "Run from git (`nextflow run lipovseklab/qc-pipeline -r <tag>`), or add --allow_dirty true for a test run."
    if( code.changed )
        error "The pipeline code has uncommitted changes:\n  ${code.changed.join('\n  ')}\n" +
              "Commit them (and tag the commit if the results may be used), or add --allow_dirty true for a test run."
    return code
}

// sha256 of a container image; cached next to the image as <image>.sha256 (keyed on size and modification time)
def imageSha256(path) {
    def f = path.toFile()
    def cache = new File("${f}.sha256")
    def key = "${f.length()} ${f.lastModified()}"
    try {
        if( cache.exists() ) {
            def lines = cache.text.readLines()
            if( lines.size() == 2 && lines[0] == key ) return lines[1]
        }
    } catch( Exception e ) {
        log.debug "Can't read ${cache.name}, hashing the image instead: ${e.message}"
    }
    def proc = ['sha256sum', f.toString()].execute()
    def out = new StringBuilder()
    def err = new StringBuilder()
    proc.waitForProcessOutput(out, err)
    if( proc.exitValue() != 0 ) return null
    def sha = out.toString().tokenize(' ')[0]
    try { cache.text = "${key}\n${sha}\n" } catch( Exception e ) { log.debug "Can't cache ${f.name} hash: ${e.message}" }
    return sha
}

// make params JSON-friendly (paths, GStrings and other objects -> text)
def plain(v) {
    if( v instanceof Map )        return v.collectEntries { k, x -> [(k.toString()): plain(x)] }
    if( v instanceof Collection ) return v.collect { x -> plain(x) }
    if( v == null || v instanceof Number || v instanceof Boolean ) return v
    return v.toString()
}

def writeProvenance(Map code, Map dataset, Map extra = [:]) {
    def info = [
        pipeline   : [name: workflow.manifest.name, version: workflow.manifest.version] + code,
        run        : [
            name        : workflow.runName,
            session_id  : workflow.sessionId?.toString(),
            started     : workflow.start?.toString(),
            command_line: workflow.commandLine,
            profile     : workflow.profile,
            launch_dir  : workflow.launchDir?.toString(),
            work_dir    : workflow.workDir?.toString(),
            user        : workflow.userName,
            stub_run    : workflow.stubRun,
            resumed     : workflow.resume,
            allow_dirty : params.allow_dirty,
        ],
        nextflow   : [version: nextflow.version.toString(), build: nextflow.build],
        dataset    : plain(dataset),
        params     : plain(params.findAll { k, _v -> k != 'datasets' }),
    ] + plain(extra)
    def out = file("${params.outdir}/pipeline_info/provenance.json")
    out.parent.mkdirs()
    out.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info)) + '\n'
    return out
}

// container images in params.container_dir with their sha256 (skipped for stub runs)
def containerHashes() {
    if( workflow.stubRun ) return [:]
    // the hashes are extra information: never let them stop the final provenance.json being written
    def dir = null
    try {
        dir = file(params.container_dir)
        if( !dir.exists() ) return [:]
        return dir.listFiles().findAll { f -> f.name.endsWith('.sif') }.sort { f -> f.name }.collectEntries { f -> [(f.name): imageSha256(f)] }
    } catch( Exception e ) {
        log.warn "Couldn't hash the container images in ${dir ?: params.container_dir}: ${e.message}"
        return [error: e.message]
    }
}
