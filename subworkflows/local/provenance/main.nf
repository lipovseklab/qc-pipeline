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

// commit, tag and uncommitted (tracked) files of the pipeline code
def codeState() {
    def commit = git(['rev-parse', 'HEAD'])
    def changed = commit ? (git(['status', '--porcelain', '--untracked-files=no']) ?: '').readLines()*.trim().findAll() : []
    [
        commit   : commit ?: workflow.commitId,
        revision : workflow.revision ?: git(['rev-parse', '--abbrev-ref', 'HEAD']),
        describe : commit ? git(['describe', '--tags', '--always', '--dirty']) : null,
        remote   : workflow.repository ?: git(['config', '--get', 'remote.origin.url']),
        changed  : changed,
    ]
}

def checkCodeState() {
    def code = codeState()
    if( workflow.stubRun || params.allow_dirty ) {
        if( params.allow_dirty && !workflow.stubRun && (code.changed || !code.commit) )
            log.warn "Running with --allow_dirty: these results can't be tied to a commit (edited files are listed in provenance.json)"
        return code
    }
    if( !code.commit )
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
    if( cache.exists() ) {
        def lines = cache.text.readLines()
        if( lines.size() == 2 && lines[0] == key ) return lines[1]
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
    def dir = file(params.container_dir)
    if( !dir.exists() ) return [:]
    dir.listFiles().findAll { f -> f.name.endsWith('.sif') }.sort { f -> f.name }.collectEntries { f -> [(f.name): imageSha256(f)] }
}
