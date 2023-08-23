#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKFORMAT } from './../../../../modules/local/checkformat'

// Test with paired-end data

workflow test_checkformat_fastagz {
    input = file(params.genomeskim_test_data['dummy']['dummy_a_fasta_gz'], checkIfExists: true)

    CHECKFORMAT(input)
}

workflow test_checkformat_fastqgz {
    input = file(params.genomeskim_test_data['dummy']['dummy_reads_1_fastq_gz'], checkIfExists: true)

    CHECKFORMAT(input)
}

workflow test_checkformat_gb {
    input = file(params.genomeskim_test_data['dummy']['dummy_genbank'], checkIfExists: true)

    CHECKFORMAT(input)
}
