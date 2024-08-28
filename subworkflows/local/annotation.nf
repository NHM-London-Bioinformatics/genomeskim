//
// Annotate contigs
//

include { MITOS                     } from '../../modules/local/mitos'
include { UNTAR as UNTAR_MITOSREFDB } from '../../modules/nf-core/untar'

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
        UNTAR_MITOSREFDB(ch_mitos_ref)

        MITOS(
            ch_contigs,
            UNTAR_MITOSREFDB.out.untar
        )
        ch_annotation_versions = ch_annotation_versions.mix(MITOS.out.versions)

        // Convert the annotations somehow?

    emit:
        //annotations = //Some annotations channel
        versions = ch_annotation_versions
}

