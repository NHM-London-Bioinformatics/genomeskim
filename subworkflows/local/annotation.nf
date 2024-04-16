//
// Annotate contigs
//

include { MITOS   } from '../modules/local/mitos'
include { BARRNAP } from '../modules/local/barrnap'

workflow ANNOTATION {
    take:
        ch_contigs
        ch_mitos_ref
        params

    main:
        ch_annotation_versions = Channel.empty()

        //
        // MODULE: MITOS
        //
        ch_mitos_ref_untar = UNTAR(ch_mitos_ref).collect()
        MITOS(
            ch_contigs,
            ch_mitos_ref_untar
        )
        ch_annotation_versions = ch_validation_versions.mix(MITOS.out.versions)

        //
        // MODULE: BARRNAP
        //
        BARRNAP(
            ch_contigs
        )
        ch_annotation_versions = ch_validation_versions.mix(BARRNAP.out.versions)


    emit:
        annotations = //Some annotations channel
        versions = ch_annotation_versions
}






}
