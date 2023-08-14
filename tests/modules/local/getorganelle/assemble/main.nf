#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GETORGANELLE } from '../../../../../modules/local/getorganelle/assemble'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_getorganelle_paired_end_animal_mt {
    input = [
        [id: 'test', single_end: false],
        [
            file("$test_data_base/mus_musculus/fastp_filtered/input_1.fastp.fastq.gz", checkIfExists: true),
            file("$test_data_base/mus_musculus/fastp_filtered/input_2.fastp.fastq.gz", checkIfExists: true),
        ]
    ]
    seeds = file("$test_data_base/get_organelle_animal_mt_0.0.1/seeds.fasta.gz", checkIfExists: true)
    labels = file("$test_data_base/get_organelle_animal_mt_0.0.1/labels.fasta.gz", checkIfExists: true)

    GETORGANELLE(input, seeds, labels)

}

