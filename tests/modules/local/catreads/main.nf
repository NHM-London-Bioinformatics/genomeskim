#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CATREADS } from '../../../../../modules/local/utilities/catreads'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_catreads_fastq {
    input = [
        [id: 'test', single_end: false],
        [
            file("$test_data_base/dummy/reads_1.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/reads_2.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/mapped_1.paired.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/mapped_2.paired.fastq.gz", checkIfExists: true)
        ]
    ]

    CATREADS(input)

}
