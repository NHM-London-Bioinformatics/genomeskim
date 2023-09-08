//
// Annotate contigs
//

include {    } from '../modules//'

workflow ANNOTATION {
    take:
        ch_contigs
        params

    main:
        ch_annotation_versions = Channel.empty()


        // If using mitos
        MITOS(ch_contigs, params.mitos_refseqdb)

        // If using Barrnap
        BARRNAP(ch_contigs)






}
