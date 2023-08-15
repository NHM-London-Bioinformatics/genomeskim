#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { SPLITREADS } from '../../../../../modules/local/getorganelle/splitreads'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_splitreads {
   input = [
        [id: 'test', single_end: false],
        [
            file("$test_data_base/dummy/reads_1.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/reads_2.fastq.gz", checkIfExists: true)
        ],
        [
            file("$test_data_base/dummy/mapped_1.paired.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/mapped_2.paired.fastq.gz", checkIfExists: true)
        ],
        [
            file("$test_data_base/dummy/mapped_1.unpaired.fastq.gz", checkIfExists: true),
            file("$test_data_base/dummy/mapped_2.unpaired.fastq.gz", checkIfExists: true)
        ]
    ]

    SPLITREADS(input)

}

