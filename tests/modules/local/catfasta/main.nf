#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CATFASTA } from '../../../../../modules/local/utilities/catfasta'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_catfasta_gz {
    input = [
        file("$test_data_base/dummy/a.fasta.gz", checkIfExists: true),
        file("$test_data_base/dummy/b.fasta.gz", checkIfExists: true)
    ]

    CATFASTA(input)
}
