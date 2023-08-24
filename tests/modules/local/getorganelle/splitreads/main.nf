#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { SPLITREADS } from './../../../../../modules/local/getorganelle/splitreads'

// Test with paired-end data

workflow test_splitreads {
    input = [
        [id: 'test', single_end: false],
        [
            file(params.genomeskim_test_data['dummy']['dummy_reads_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_reads_2_fastq_gz'], checkIfExists: true)
        ],
        [
            file(params.genomeskim_test_data['dummy']['dummy_mapped_paired_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_mapped_paired_2_fastq_gz'], checkIfExists: true)
        ],
        [
            file(params.genomeskim_test_data['dummy']['dummy_mapped_unpaired_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['dummy']['dummy_mapped_unpaired_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    SPLITREADS(input)

}

