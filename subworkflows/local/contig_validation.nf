//
// Sort out all the different possible reference routes and do any combining needed
//

include { GOFETCH                     } from '../modules/local/gofetch'

workflow CONTIG_VALIDATION {
    take:
        params

    main:
        //
        // MODULE: BLASTN
        //
        BLAST_BLASTNREMOTE(GETORGANELLE.out.contig)

        //
        // MODULE: MINIMAP2
        //
        MINIMAP2_ALIGN()

    // OUTPUT
    emit:

}
