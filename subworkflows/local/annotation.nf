//
// Annotate contigs
//

include { MITOS   } from '../../modules/local/mitos'
include { BARRNAP } from '../../modules/local/barrnap'
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

        //
        // MODULE: BARRNAP
        //
/*         BARRNAP(
            ch_contigs
        )
        ch_annotation_versions = ch_annotation_versions.mix(BARRNAP.out.versions)
*/

    emit:
        annotations = //Some annotations channel
        versions = ch_annotation_versions
}

