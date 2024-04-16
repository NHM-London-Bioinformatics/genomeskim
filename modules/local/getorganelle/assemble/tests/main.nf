#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GETORGANELLE } from './../../../../../modules/local/getorganelle/assemble'

// Test with paired-end data

workflow test_getorganelle_paired_end_animal_mt {
    input = [
        [id: 'test', single_end: false],
        [
            file(params.genomeskim_test_data['mus_musculus']['fastp_filtered_1_fastq_gz'], checkIfExists: true),
            file(params.genomeskim_test_data['mus_musculus']['fastp_filtered_2_fastq_gz'], checkIfExists: true),
        ]
    ]
    seeds = file(params.genomeskim_test_data['references']['getorganelle_animal_mt_seeds_fasta_gz'], checkIfExists: true)
    labels = file(params.genomeskim_test_data['references']['getorganelle_animal_mt_labels_fasta_gz'], checkIfExists: true)

    GETORGANELLE(input, seeds, labels)

}

