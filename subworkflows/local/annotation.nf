//
// Annotate contigs
//

include { MITOS   } from '../../modules/local/mitos'
include { UNTAR   } from '../../modules/nf-core/untar'

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
        UNTAR(ch_mitos_ref)

        MITOS(
            ch_contigs,
            UNTAR.out.untar
        )
        ch_annotation_versions = ch_annotation_versions.mix(MITOS.out.versions)

        // Convert the annotations somehow?

    emit:
        //annotations = //Some annotations channel
        versions = ch_annotation_versions
}

