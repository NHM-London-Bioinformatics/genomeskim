#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GETGOREFS } from '../../../../../modules/local/getorganelle/getgoreferences'

// Test with paired-end data

test_data_base = '../../../../../test_data'

workflow test_getgoreferences_animal_mt {

    GETGOREFS('animal_mt')

}

