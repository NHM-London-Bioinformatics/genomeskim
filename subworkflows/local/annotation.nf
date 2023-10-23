//
// Annotate contigs
//

include {    } from '../modules//'

workflow ANNOTATION {
    take:
        ch_contigs
        ch_mitos_ref
        params

    main:
        ch_annotation_versions = Channel.empty()

        // If using mitos
        ch_mitos_ref_untar = UNTAR(ch_mitos_ref).collect()
        MITOS(ch_contigs, ch_mitos_ref_untar)

        // If using Barrnap
        BARRNAP(ch_contigs)






}
