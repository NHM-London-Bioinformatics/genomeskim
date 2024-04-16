#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CATREADS } from './../../../../modules/local/catreads'

// Test with paired-end data

workflow test_catreads_fastq {
    input = [
        [id: 'test', single_end: false],
        [
            file(params.genomeskim_test_data['dummy']['dummy_reads_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_reads_2_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_mapped_paired_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_mapped_paired_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    CATREADS(input, 'test')

}
