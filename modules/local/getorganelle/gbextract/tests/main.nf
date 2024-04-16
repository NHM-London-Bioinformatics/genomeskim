#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GBEXTRACT } from './../../../../../modules/local/getorganelle/gbextract'

// Test with paired-end data

workflow test_gbextract {

    input = file(params.genomeskim_test_data['dummy']['dummy_genbank'], checkIfExists: true)

    GBEXTRACT(input)

}

