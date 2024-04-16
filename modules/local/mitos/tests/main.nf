#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOS } from './../../../../modules/local/mitos'

// Test with paired-end data

workflow test_mitos_m5 {
    contig = [
        [id: 'test'],
        [file(params.genomeskim_test_data['dummy']['dummy_annotate_m5_fasta_gz'], checkIfExists: true)]
    ]
    db = file(params.mitos_ref_databases['refseq89m']['file'])

    MITOS(contig, db)
}
