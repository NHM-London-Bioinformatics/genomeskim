//
// Validate and filter organelle contigs
//

include { BLAST_BLASTN           } from '../../modules/nf-core/blast/blastn'
include { MINIMAP2_ALIGN         } from '../../modules/nf-core/minimap2/align'
include { BLOBTOOLS              } from '../../modules/local/blobtools'
include { UNTAR as UNTAR_TAXDUMP } from '../../modules/nf-core/untar'
include { UNTAR as UNTAR_BLASTDB } from '../../modules/nf-core/untar'

workflow VALIDATION {
    take:
        ch_reads
        ch_contigsin
        ch_taxdump
        ch_blastdbin
        params

    main:
        ch_validation_versions = Channel.empty()

        // Prepare taxdump and blastdb
        UNTAR_TAXDUMP(ch_taxdump)

        // Prepare BLAST DB
        if ( params.blastdbpath ){
            ch_blastdb = ch_blastdbin
        } else {
            UNTAR_BLASTDB(ch_blastdbin).out.untar.set(ch_blastdb)
        }


        // Split out contigs to multiple parallel channels

        ch_contigssplit = ch_contigsin
            .multiMap { i -> contigs4blast: contigs4map: contigs4blob: i }

        //
        // MODULE: BLASTN
        //
        BLAST_BLASTN(
            ch_contigssplit.contigs4blast,
            ch_blastdb
        )
        ch_validation_versions = ch_validation_versions.mix(BLAST_BLASTN.out.versions)

        //
        // MODULE: MINIMAP2
        //

        // Join contigs and reads to sync order, then split to form parallel channels
        ch_mapinput = ch_contigssplit.contigs4map
            .join(ch_reads)
            .multiMap{ i ->
                contigs: [i[0], i[1]]
                reads:   [i[0], i[2]]
            }

        MINIMAP2_ALIGN(
            ch_mapinput.reads,
            ch_mapinput.contigs,
            true, // BAM format output
            false, // BAM index extension
            false, // CIGAR PAF format
            false  // CIGAR BAM
        )
        ch_validation_versions = ch_validation_versions.mix(MINIMAP2_ALIGN.out.versions)

        //
        // MODULE: BLOBTOOLS
        //

        // Join all results to sync order, then split to form parallel channels
        ch_blobinput = ch_contigssplit.contigs4blob
            .join(BLAST_BLASTN.out.txt)
            .join(MINIMAP2_ALIGN.out.bam)
            .multiMap{ i ->
                contigs: [i[0], i[1]]
                blast:   [i[0], i[2]]
                mapd:    [i[0], i[3]]
            }

        BLOBTOOLS(
            ch_blobinput.contigs,
            ch_blobinput.blast,
            ch_blobinput.mapd,
            UNTAR_TAXDUMP.out.untar
        )
        ch_validation_versions = ch_validation_versions.mix(BLOBTOOLS.out.versions)


    // OUTPUT
    emit:
        versions = ch_validation_versions

}
