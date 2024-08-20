process GETORGANELLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.1--pyhdfd78af_0':
        'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(reads)
        path seeds
        path genes

    output:
        //TODO: MULTIQC doesn't support getorganelle
        tuple val(meta), path("output/*1.1.path_sequence.fasta.gz")  , emit: contigs
        tuple val(meta), path("output/extended_*_paired.fq.gz")      , emit: pairedreads
        tuple val(meta), path("output/extended_*_unpaired.fq.gz")    , emit: unpairedreads
        tuple val(meta), path("getorganelle.log")                    , emit: log
        path "versions.yml"                                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """

        zcat $seeds > seeds.fasta
        zcat $genes > labels.fasta

        get_organelle_from_reads.py \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o output \\
            -t $task.cpus \\
            --verbose \\
            -s seeds.fasta \\
            --genes labels.fasta \\
            $args \\
            2>&1 > getorganelle.log

        rm seeds.fasta labels.fasta

        gzip output/*.fasta
        gzip output/*.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            getorganelle: \$(get_organelle_from_reads.py -v 2>&1 | sed -e "s/^.* v//g")
        END_VERSIONS

        """

    stub:
        def args = task.ext.args ?: ''
        """
        mkdir output

        touch output/contigs1.1.path_sequence.fasta.gz
        touch output/extended_{1,2}_paired.fq.gz
        touch output/extended_{1,2}_unpaired.fq.gz
        touch getorganelle.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            getorganelle: \$(get_organelle_from_reads.py -v 2>&1 | sed -e "s/^.* v//g")
        END_VERSIONS
        """

}

