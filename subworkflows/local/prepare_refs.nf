//
// Sort out all the different possible reference routes and do any combining needed
//

include { GOFETCH                     } from '../modules/local/gofetch'
include { CHECKFORMAT                 } from '../modules/local/utilities/checkformat'
include { GBEXTRACT                   } from '../modules/local/getorganelle/gbextract'
include { GETGOREFS                   } from '../modules/local/getorganelle/getgoreferences'
include { CATFASTAORG                 } from '../modules/local/utilities/catfasta'
include { CATFASTAGENE                } from '../modules/local/utilities/catfasta'
include { CATFASTASEED                } from '../modules/local/utilities/catfasta'
include { CATFASTALABELS              } from '../modules/local/utilities/catfasta'

workflow PREPARE_REFS {
    take:
        params

    main:
        ch_preprefs_versions = Channel.empty()

        ch_orgseqs = Channel.empty()
        ch_orggenes = Channel.empty()

        ch_goseeds = Channel.empty()
        ch_golabels = Channel.empty()

        ch_mitos_ref = Channel.empty()

        // If taxon or lineage are supplied
        if ( params.gofetch_taxon || params.gofetch_lineage ){
            GOFETCH(params.gofetch_taxon, params.gofetch_lineage)

            ch_orgseqs = ch_orgseqs.mix(GOFETCH.out.seqs)
            ch_orggenes = ch_orggenes.mix(GOFETCH.out.genes)
            ch_preprefs_versions = ch_preprefs_versions.mix(GOFETCH.out.versions)

        }

        // If reference is supplied
        if ( params.organellerefs ){
            ch_refpath = Channel.fromPath(params.organellerefs).collect()

            CHECKFORMAT(ch_refpath)

            if( CHECKFORMAT.out.format == 'gb' ){

                GBEXTRACT(ch_refpath)

                ch_orgseqs = ch_orgseqs.mix(GBEXTRACT.out.seqs)
                ch_orggenes = ch_orggenes.mix(GBEXTRACT.out.genes)
                ch_preprefs_versions = ch_preprefs_versions.mix(GBEXTRACT.out.versions)

            } elif ( CHECKFORMAT.out.format == 'fasta' ){

                ch_orgseqs = ch_orgseqs.mix(ch_refpath)

            } else { exit 1, 'Organelle references file is not fasta or genbank flat file format'}

        }

        // Merge files
        ch_orgseqs  = CATFASTAORG(ch_orgseqs.collect()).out.path.collect()
        ch_orggenes = CATFASTAGENE(ch_orggenes.collect()).out.path.collect()

        // Get GetOrganelle references if needed, and then merge as necessary
        if ( params.getorganelle_ref_action ) {

            GETGOREFS(params.getorganelle_genometype)
            ch_preprefs_versions = ch_preprefs_versions.mix(GETGOREFS.out.versions)

            // Store the seeds and/or labels in the relevant channels if they're being used
            // params.getorganelle_ref_action should be only_both if non of taxon/lineage/reference are supplied
            if ( params.getorganelle_ref_action in ['add_seeds', 'add_both', 'only_seeds', 'only_both'] ) {
                ch_goseeds = ch_goseeds.mix(GETGOREFS.out.seeds)
            }
            if ( params.getorganelle_ref_action in ['add_labels', 'add_both', 'only labels', 'only_both'] ) {
                ch_golabels = ch_golabels.mix(GETGOREFS.out.labels)
            }

            // Add in formatted references if needed
            if ( params.getorganelle_ref_action in ['add_seeds', 'add_both'] ) {
                ch_goseeds = ch_goseeds.mix(ch_orgseqs)
            }
            if ( params.getorganelle_ref_action in ['add_labels', 'add_both'] ) {
                ch_golabels = ch_golabels.mix(ch_orggenes)
            }

            ch_goseeds = CATFASTASEED(ch_goseeds.collect()).out.path.collect()
            ch_golabels = CATFASTALABELS(ch_golabels.collect()).out.path.collect()

        }

        // If MITOS is to be run
        if ( params.mitos_refdbid ) {

            ch_mitos_ref = Channel.fromPath(params.mitos_ref_databases[params.mitos_refdbid]["file"]).collect()

        }

    // OUTPUT
    emit:
        orgseqs  = ch_orgseqs
        goseeds  = ch_goseeds
        golabels = ch_golabels
        mitosref = ch_mitos_ref
        versions = ch_preprefs_versions
}
