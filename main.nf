nextflow.enable.dsl = 2

params.samplesheet = null
params.chm13_fasta = null
params.hor_consensuses = "${projectDir}/22q11.2_25-28-41-105-120_HORs.fasta"
params.annotation = "${projectDir}/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"
params.reference_record = 'CHM13#0#chr22:18000000-23000000'
params.outdir = 'results'

process BUILD_FAMILY_FASTA {
    tag { family }

    publishDir "${params.outdir}/family_fastas", mode: 'copy', overwrite: false,
        saveAs: { filename -> "${family}/${filename}" }

    input:
    tuple val(family), path(assemblies)
    path chm13_fasta

    output:
    tuple val(family), path('family.fa'), emit: family_fasta

    script:
    def assembly_files = assemblies.collect { "'${it}'" }.join(' ')

    """
    cat "${chm13_fasta}" ${assembly_files} > family.fa
    """
}

process WFMASH_FAMILY_ALIGNMENT {
    tag { family }
    cpus 64

    publishDir "${params.outdir}/alignments", mode: 'copy', overwrite: false,
        saveAs: { filename -> "${family}/${filename}" }

    input:
    tuple val(family), path(family_fasta)

    output:
    tuple val(family), path('alignments.paf'), emit: paf

    """
    wfmash -t 64 -p 95 -Y '#' -S 20k -j 200k --scaffold-out \\
        "${family_fasta}" > alignments.paf
    """
}

process QUERY_CHR22Q_INTERVALS {
    tag { family }

    publishDir "${params.outdir}/chr22q_intervals", mode: 'copy', overwrite: false,
        saveAs: { filename -> "${family}/${filename}" }

    input:
    tuple val(family), path(paf)

    output:
    tuple val(family), path('chr22q_intervals.txt'), emit: intervals

    """
    impg query -a "${paf}" -r "${params.reference_record}" -d 1000000 |
        awk '\$3-\$2 > 100e3 {print \$1":"\$2+1"-"\$3}' > chr22q_intervals.txt
    """
}

process EXTRACT_CHR22Q_FASTA {
    tag { family }

    publishDir "${params.outdir}/chr22q_fastas", mode: 'copy', overwrite: false,
        saveAs: { filename -> "${family}/${filename}" }

    input:
    tuple val(family), path(family_fasta), path(intervals)

    output:
    tuple val(family), path('chr22q.fa'), emit: chr22q_fasta

    """
    samtools faidx "${family_fasta}" -r "${intervals}" > chr22q.raw.fa
    awk -v reference="${params.reference_record}" \\
        'index(\$0, ">" reference ":") == 1 { print ">" reference; next } { print }' \\
        chr22q.raw.fa > chr22q.fa
    """
}

process RUN_CHR22Q_LCR_ANALYSIS {
    tag { family }
    cpus 8

    publishDir "${params.outdir}/lcr_analysis", mode: 'copy', overwrite: false,
        saveAs: { filename -> "${family}/${filename}" }

    input:
    tuple val(family), path(chr22q_fasta)
    path hor_consensuses
    path annotation

    output:
    tuple val(family), path('chr22q_lcr_analysis'), emit: analysis

    script:
    """
    "${projectDir}/run_chr22q_lcr_analysis.sh" \\
        --fasta "${chr22q_fasta}" \\
        --hor-consensuses "${hor_consensuses}" \\
        --annotation "${annotation}" \\
        --reference-record "${params.reference_record}" \\
        --output chr22q_lcr_analysis
    """
}

workflow {
    if (!params.samplesheet) {
        error 'Provide a two-column CSV sample sheet with --samplesheet.'
    }
    if (!params.chm13_fasta) {
        error 'Provide the CHM13 FASTA with --chm13_fasta.'
    }

    family_assemblies = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.ASSEMBLY || !row.FAMILY) {
                error 'The sample sheet must contain non-empty ASSEMBLY and FAMILY columns.'
            }
            tuple(row.FAMILY.toString(), file(row.ASSEMBLY.toString(), checkIfExists: true))
        }
        .groupTuple()

    chm13_fasta = Channel.value(file(params.chm13_fasta, checkIfExists: true))

    family_fastas = BUILD_FAMILY_FASTA(family_assemblies, chm13_fasta)
    family_alignments = WFMASH_FAMILY_ALIGNMENT(family_fastas)
    family_intervals = QUERY_CHR22Q_INTERVALS(family_alignments)
    chr22q_inputs = family_fastas.join(family_intervals)
    chr22q_fastas = EXTRACT_CHR22Q_FASTA(chr22q_inputs)

    hor_consensuses = Channel.value(file(params.hor_consensuses, checkIfExists: true))
    annotation = Channel.value(file(params.annotation, checkIfExists: true))
    RUN_CHR22Q_LCR_ANALYSIS(chr22q_fastas, hor_consensuses, annotation)
}
