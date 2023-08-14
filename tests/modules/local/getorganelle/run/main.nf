#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GETORGANELLE } from '../../../../../modules/local/getorganelle/run'

// Test with paired-end data

test_data_base = '../../../../../test_data/modules/local/getorganelle/run'

workflow test_getorganelle_paired_end {
    input = [
        [id: 'test', single_end: false],
        [
            file("$test_data_base/input_1.fastq.gz"),
            file("$test_data_base/input_2.fastq.gz"),
        ]
    ]
    seeds = file("$test_data_base/seeds.fasta.gz")
    labels = file("$test_data_base/labels.fasta.gz")

    GETORGANELLE(input, seeds, labels)

}

