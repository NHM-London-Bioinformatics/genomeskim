#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CATFASTA } from './../../../../modules/local/catfasta/main.nf'

// Test with paired-end data

workflow test_catfasta_gz {
    input = [
        file(params.genomeskim_test_data['dummy']['dummy_a_fasta_gz'], checkIfExists: true),
        file(params.genomeskim_test_data['dummy']['dummy_b_fasta_gz'], checkIfExists: true)
    ]

    CATFASTA(input)
}
