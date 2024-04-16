//
// Validate and filter organelle contigs
//

include { BLAST_BLASTN     } from '../../modules/nf-core/blast/blastn'
include { MINIMAP2_ALIGN   } from '../../modules/nf-core/minimap2/align'
include { BLOBTOOLS        } from '../../modules/local/blobtools'

workflow ORGANELLE_VALIDATION {
    take:
        ch_reads
        ch_contigsin
        params

    main:
        ch_validation_versions = Channel.empty()


        // Split out contigs to multiple parallel channels

        ch_contigssplit = ch_contigsin
            .multiMap { i -> contigs4blast: contigs4map: contigs4blob: i}

        //
        // MODULE: BLASTN
        //
        BLAST_BLASTN(
            ch_contigssplit.contigs4blast,
            params.blastdbpath
        )
        ch_validation_versions = ch_validation_versions.mix(BLAST_BLASTN.out.versions)

        //
        // MODULE: MINIMAP2
        //
        MINIMAP2_ALIGN(
            ch_reads
            ch_contigssplit.contigs4map,
            true, // BAM format output
            false, // CIGAR PAF format
            false  // CIGAR BAM
        )
        ch_validation_versions = ch_validation_versions.mix(MINIMAP2_ALIGN.out.versions)

        // MODULE: BLOBTOOLS
        BLOBTOOLS(
            ch_contigs4blob
            BLAST_BLASTN.out.txt
            MINIMAP2_ALIGN.out.bam
            params.taxdump
        )
        ch_validation_versions = ch_validation_versions.mix(BLOBTOOLS.out.versions)

        // TODO: Some filtering somewhere to remove contigs that we don't want

    // OUTPUT
    emit:
        contigs  = //Some contigs channel
        versions = ch_validation_versions
}
