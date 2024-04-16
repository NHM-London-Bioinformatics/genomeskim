/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def addtomap(map, key, value) {
    nwmap = map.clone()
    nwmap.put(key, value)
    return nwmap
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GETORGANELLE                } from '../modules/local/getorganelle/assemble'
include { SPLITREADS                  } from '../modules/local/getorganelle/splitreads/main'
include { CATREADS                    } from '../modules/local/catreads'
include { GENOMESCOPE2                } from '../modules/local/genomescope2/main'

//TODO modules for art, skmer


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// This subworkflow is all about parsing and checking the samplesheet
//include { INPUT_CHECK          } from '../subworkflows/local/input_check'
include { PREPARE_REFS         } from '../subworkflows/local/prepare_refs'
include { ORGANELLE_VALIDATION } from '../subworkflows/local/organelle_validation'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Installed from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomeskim_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMESKIM {

    take:
    ch_samplesheet

    main:
    // Create an empty channel to record software versions
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    /* INPUT_CHECK (
        file(params.input)
    ) */
    // mix adds the contents of another channel (in this case out.versions from INPUT_CHECK) to
    // the given channel (ch_versions)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // SUBWORKFLOW: Retrieve reference sequences as determined by the user
    //

    PREPARE_REFS(params)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

//
    // MODULE: Run fastp
    //
    FASTP (
        INPUT_CHECK.out.reads,
        adapters, // Just passing a value, no need for this to be a channel
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_cleanreads = FASTP.out.reads.
        multiMap{ it -> goreads: binreads: mapreads: it }

    //
    // MODULE: GetOrganelle
    //

    GETORGANELLE (
        ch_cleanreads.goreads,
        PREPARE_REFS.out.goseeds,
        PREPARE_REFS.out.golabels
    )
    // Split the input reads based on the files comprising the paired reads and unpaired reads used by getorganelle
    ch_getorganelle_readuse = FASTP.out.reads.join(GETORGANELLE.out.pairedreads)
    ch_getorganelle_readuse = ch_getorganelle_readuse.join(GETORGANELLE.out.unpairedreads)

    SPLITREADS (
        ch_getorganelle_readuse
    )

    ch_versions = ch_versions.mix(GETORGANELLE.out.versions.first())
    ch_versions = ch_versions.mix(SPLITREADS.out.versions.first())

    // Concatenate unpaired and unused reads
    ch_nucreads = SPLITREADS.out.usedreadsup
        .mix(SPLITREADS.out.unusedreads)
        .groupTuple().map { i -> [ i[0], i[1].flatten() ] }

    CATREADS (
        ch_nucreads,
        'nuclear'
    )

    //
    // Run contig validation
    //
    ORGANELLE_VALIDATION(
        ch_cleanreads.mapreads,
        GETORGANELLE.out.contigs,
        params
    )
    ch_versions = ch_versions.mix(ORGANELLE_VALIDATION.out.versions)

    //
    // Run annotation
    //
    ANNOTATION(
        ORGANELLE_VALIDATION.out.contigs,
        ch_mitos_ref,
        params
    )

    //
    // Extract barcodes
    //


    //
    // MODULE: GENOMESCOPE2
    //
    GENOMESCOPE2(CATNUCREADS.out.catreads)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
