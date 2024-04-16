#!/usr/bin/env nextflow

include { GETGOREFS } from './../../../../../modules/local/getorganelle/getgoreferences'

// Test with paired-end data

workflow test_getgoreferences_animal_mt {

    GETGOREFS('animal_mt')

}

