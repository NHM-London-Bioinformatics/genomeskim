//
// Validate and filter organelle contigs
//

include { BLAST_BLASTN           } from '../../modules/nf-core/blast/blastn'
include { MINIMAP2_ALIGN         } from '../../modules/nf-core/minimap2/align'
include { BLOBTOOLS              } from '../../modules/local/blobtools'
include { UNTAR as UNTAR_TAXDUMP } from '../../modules/nf-core/untar'

workflow ORGANELLE_VALIDATION {
    take:
        ch_reads
        ch_contigsin
        ch_blastdbpath
        ch_taxdump
        params

    main:
        ch_validation_versions = Channel.empty()


        // Split out contigs to multiple parallel channels

        ch_contigssplit = ch_contigsin
            .multiMap { i -> contigs4blast: contigs4map: contigs4blob: i }

        //
        // MODULE: BLASTN
        //
        BLAST_BLASTN(
            ch_contigssplit.contigs4blast,
            ch_blastdbpath
        )
        ch_validation_versions = ch_validation_versions.mix(BLAST_BLASTN.out.versions)

        //
        // MODULE: MINIMAP2
        //
        MINIMAP2_ALIGN(
            ch_reads,
            ch_contigssplit.contigs4map,
            true, // BAM format output
            false, // BAM index extension
            false, // CIGAR PAF format
            false  // CIGAR BAM
        )
        ch_validation_versions = ch_validation_versions.mix(MINIMAP2_ALIGN.out.versions)

        //
        // MODULE: BLOBTOOLS
        //
        UNTAR_TAXDUMP(ch_taxdump)

        BLOBTOOLS(
            ch_contigssplit.contigs4blob,
            BLAST_BLASTN.out.txt,
            MINIMAP2_ALIGN.out.bam,
            UNTAR_TAXDUMP.out.untar
        )
        ch_validation_versions = ch_validation_versions.mix(BLOBTOOLS.out.versions)

    // OUTPUT
    emit:
        versions = ch_validation_versions
}
