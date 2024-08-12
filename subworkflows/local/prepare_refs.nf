//
// Sort out all the different possible reference routes and do any combining needed
//

include { GOFETCH                    } from '../../modules/local/gofetch'
include { CHECKFORMAT                } from '../../modules/local/checkformat'
include { GBEXTRACT                  } from '../../modules/local/getorganelle/gbextract'
include { GETGOREFS                  } from '../../modules/local/getorganelle/getgoreferences'
include { CATFASTA as CATFASTAORG    } from '../../modules/local/catfasta'
include { CATFASTA as CATFASTAGENE   } from '../../modules/local/catfasta'
include { CATFASTA as CATFASTASEED   } from '../../modules/local/catfasta'
include { CATFASTA as CATFASTALABEL  } from '../../modules/local/catfasta'

workflow PREPARE_REFS {
    take:
        params

    main:
        ch_preprefs_versions = Channel.empty()

        ch_orgseqs = Channel.empty()
        ch_orggenes = Channel.empty()

        ch_goseeds = Channel.empty()
        ch_golabels = Channel.empty()

        referenceformat = ''
        nocustom = !params.gofetch_taxon & !params.gofetch_lineage & !params.organellerefs

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

            referenceformat = CHECKFORMAT(ch_refpath).out.format

            if( referenceformat == 'gb' ){

                GBEXTRACT(ch_refpath)

                ch_orgseqs = ch_orgseqs.mix(GBEXTRACT.out.seqs)
                ch_orggenes = ch_orggenes.mix(GBEXTRACT.out.genes)
                ch_preprefs_versions = ch_preprefs_versions.mix(GBEXTRACT.out.versions)

            } elif ( referenceformat == 'fasta' ){

                ch_orgseqs = ch_orgseqs.mix(ch_refpath)

            } else { error('Organelle references file is not fasta or genbank flat file format')}

        }

        // Merge files
        if ( params.gofetch_taxon || params.gofetch_lineage || params.organellerefs ){

            CATFASTAORG(ch_orgseqs.collect(), "org")
            CATFASTAGENE(ch_orggenes.collect(), "gene")
            ch_orgseqs  = CATFASTAORG.out.catfasta.collect()
            ch_orggenes = CATFASTAGENE.out.catfasta.collect()
        }

        // Get GetOrganelle references if needed, and then merge as necessary
        if ( params.getorganelle_ref_action || nocustom ) {
            GETGOREFS(
                Channel.fromList(params.getorganelle_genometype.tokenize(','))
            )

            ch_preprefs_versions = ch_preprefs_versions.mix(GETGOREFS.out.versions)

            // Store the seeds and/or labels in the relevant channels if they're being used
            // params.getorganelle_ref_action should be only_both if none of taxon/lineage/reference are supplied
            if ( params.getorganelle_ref_action in ['add_seeds', 'add_both', 'only_seeds', 'only_both'] || nocustom ) {

                ch_goseeds = ch_goseeds.mix(GETGOREFS.out.seeds)
            }
            if ( params.getorganelle_ref_action in ['add_labels', 'add_both', 'only_labels', 'only_both'] || nocustom ) {
                ch_golabels = ch_golabels.mix(GETGOREFS.out.labels)

            }

            // Add in formatted references if needed
            if ( params.getorganelle_ref_action in ['add_seeds', 'add_both'] ) {
                ch_goseeds = ch_goseeds.mix(ch_orgseqs)
            }
            if ( params.getorganelle_ref_action in ['add_labels', 'add_both'] ) {
                ch_golabels = ch_golabels.mix(ch_orggenes)
            }


            CATFASTASEED(ch_goseeds.collect(), "seed")
            CATFASTALABEL(ch_golabels.collect(), "label")
            ch_goseeds = CATFASTASEED.out.catfasta.collect()
            ch_golabels = CATFASTALABEL.out.catfasta.collect()


        } else {
            ch_goseeds = ch_orgseqs
            if ( params.gofetch_taxon || params.gofetch_lineage || referenceformat == 'gb' ){
                ch_golabels = ch_orggenes
            } else {
                def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    " Error: the parameter combination supplied has not resulted in any gene  \n" +
                    " sequences. Either supply a taxon or lineage to retrieve from GoFetch,   \n" +
                    " specify reference(s) in genbank format, or set --getorganelle_ref_action\n" +
                    " to add_labels, add_both or only_both to use the GetOrganelle defaults   \n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                error(error_string)
            }

        }

    // OUTPUT
    emit:
        orgseqs  = ch_orgseqs
        goseeds  = ch_goseeds
        golabels = ch_golabels
        versions = ch_preprefs_versions
}
