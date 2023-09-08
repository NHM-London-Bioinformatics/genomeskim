//
// Sort out all the different possible reference routes and do any combining needed
//

include { GOFETCH                     } from '../modules/local/gofetch'

workflow CONTIG_VALIDATION {
    take:
        params
        contig

    main:
        //
        // MODULE: BLASTN
        //
        BLAST_BLASTN(contig, params.blastdb)

        //
        // MODULE: MINIMAP2
        //
        MINIMAP2_ALIGN()

    // OUTPUT
    emit:

}
