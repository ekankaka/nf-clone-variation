nextflow.enable.dsl = 2

process ASSESS_CLONE_VARIATION {
    tag "$header"
    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/plots", mode: 'copy', pattern: '*.png', overwrite: true
    publishDir "${params.outdir}/plots", mode: 'copy', pattern: '*.pdf', overwrite: true

    input:
    tuple val(header), val(row), val(assignments), val(missing), path(alignment, stageAs: 'alignment.fasta')
    path analysis_script

    output:
    path 'report', emit: report
    path '*.png', optional: true
    path '*.pdf', optional: true

    script:
    def metadata = groovy.json.JsonOutput.toJson([clone: row, assignments: assignments, missing_alignment: missing])
    """
    cat > metadata.json <<'NF_METADATA'
    $metadata
    NF_METADATA
    python '${analysis_script}' analyze --metadata metadata.json --alignment alignment.fasta --family-fraction ${params.min_family_fraction}
    """
}

process MERGE_RESULTS {
    cpus 1
    memory '2 GB'

    publishDir params.outdir, mode: 'copy', pattern: 'clone_variation_summary.csv', overwrite: true
    publishDir params.outdir, mode: 'copy', pattern: 'variants*.csv', overwrite: true
    publishDir "${params.outdir}/review", mode: 'copy', pattern: 'skipped_clones.csv', overwrite: true

    input:
    path manifest
    path(reports, stageAs: 'clone??????/*')
    path analysis_script

    output:
    path 'clone_variation_summary.csv'
    path 'variants.csv'
    path 'variants_by_timepoint.csv'
    path 'skipped_clones.csv'

    script:
    """
    python '${analysis_script}' merge --input '${manifest}'
    """
}

workflow {
    if (!params.input || !params.clone_alignment_folder || !params.read_assignments)
        error 'Required: --input, --clone_alignment_folder and --read_assignments'
    if (params.min_unique_reads.toString().toInteger() < 3)
        error '--min_unique_reads must be at least 3'
    if (params.min_family_fraction.toString().toDouble() <= 0.5 ||
        params.min_family_fraction.toString().toDouble() > 1)
        error '--min_family_fraction must be > 0.5 and <= 1'

    def manifest = file(params.input, checkIfExists: true)
    def assignments_file = file(params.read_assignments, checkIfExists: true)
    def folder = file(params.clone_alignment_folder, checkIfExists: true)
    def analysis_script = file("$projectDir/bin/assess_clone_variation.py", checkIfExists: true)
    if (!folder.isDirectory()) error '--clone_alignment_folder must be a directory'

    // Group metadata by patient and clone, not by nucleotide sequence.
    clones = Channel.fromPath(manifest)
        .splitCsv(header: true)
        .filter { row ->
            row.n_unique_reads.toInteger() >= params.min_unique_reads.toString().toInteger()
        }
        .map { row -> tuple([row.patient_id, row.clone_id], row) }

    assignments = Channel.fromPath(assignments_file)
        .splitCsv(header: true)
        .filter { row -> row.classification == 'CLONE' }
        .map { row -> tuple([row.patient_id, row.clone_id], row) }
        .groupTuple()

    inputs = clones.join(assignments, remainder: true)
        .filter { key, row, reads -> row != null }
        .map { key, row, reads ->
            def header = row.consensus_header
            if (!header || !(header ==~ /[A-Za-z0-9][A-Za-z0-9._%~-]*/))
                error "Unsafe consensus_header: $header"
            def alignment = folder.resolve(header + '_alignment.fasta')
            if (!alignment.isFile()) alignment = folder.resolve(header + '.fasta')
            def missing = !alignment.isFile()
            // A placeholder allows missing alignments to be reported rather than lost.
            tuple(header, row, reads ?: [], missing, missing ? manifest : alignment)
        }

    ASSESS_CLONE_VARIATION(inputs, Channel.value(analysis_script))

    // Keep a staged placeholder so the merge also runs with no eligible clones.
    reports = ASSESS_CLONE_VARIATION.out.report.collect().ifEmpty([])
        .map { paths -> [manifest] + paths }
    MERGE_RESULTS(Channel.value(manifest), reports, Channel.value(analysis_script))
}
