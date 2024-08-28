process BLOBTOOLS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blobtoolkit:4.3.11':
        'genomehubs/blobtoolkit:4.3.11' }"

    input:
        tuple val(meta),        path(contig)
        tuple val(metablast),   path(blastn)
        tuple val(metamap),     path(bam)
        tuple val(metataxdump), path(taxdump)

    output:
        tuple val(meta), path('*_summary.tsv'), emit: table
        //TODO: MultiQC - currently not supported
        path "versions.yml"                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def filargs = task.ext.filargs ?: ''
        def creargs = task.ext.crargs ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """

        blobtools create \
            --taxdump $taxdump \
            --fasta $contig \
            --hits $blastn \
            --cov $bam \
            $filargs \
            blobout/ &> ${prefix}_blobtools_create.log

        blobtools filter \
            $creargs \
            --table ${prefix}_summary.tsv \
            blobout &> ${prefix}_blobtools_filter.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: \$(blobtools --version | sed -e "s/blobtoolkit //")
        END_VERSIONS
        """

    stub:
        def filargs = task.ext.filargs ?: ''
        def creargs = task.ext.crargs ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir blobout

        touch bestsumorder_positions.json
        touch gc.json
        touch identifiers.json
        touch length.json
        touch meta.json
        touch ncount.json
        touch ${prefix}_cov.json
        for i in superkingdom kingdom phylum class order family genus species
        do
            for j in _cindex '' _positions _score
            do
                touch bestsumorder_\${i}\${j}.json
            done
        done

        mv *.json blobout/

        touch ${prefix}_blobtools_create.log
        touch ${prefix}_summary.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: \$(blobtools --version | sed -e "s/blobtoolkit //")
        END_VERSIONS
        """
}
