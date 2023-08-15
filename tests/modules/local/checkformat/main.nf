#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKFORMAT } from '../../../../../modules/local/utilities/checkformat'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_checkformat_fastagz {
    input = file("$test_data_base/dummy/a.fasta.gz", checkIfExists: true)

    CHECKFORMAT(input)
}

workflow test_checkformat_fastqgz {
    input = file("$test_data_base/dummy/reads_1.fastq.gz", checkIfExists: true)

    CHECKFORMAT(input)
}

workflow test_checkformat_gb {
    input = file("$test_data_base/dummy/dummy.gb", checkIfExists: true)

    CHECKFORMAT(input)
}
