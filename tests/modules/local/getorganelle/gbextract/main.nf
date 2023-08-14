#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GBEXTRACT } from '../../../../../modules/local/getorganelle/gbextract'

// Test with paired-end data

test_data_base = '../../../../../test_data/'

workflow test_gbextract {

    input = file("$test_data_base/dummy/dummy.gb", checkIfExists: true)

    GBEXTRACT(input)

}

