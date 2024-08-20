process JELLYFISH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::jellyfish==2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jellyfish:2.2.3--h6bb024c_3':
        'biocontainers/jellyfish:v2.2.3_cv3' }"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.jf"),    emit: jf
        tuple val(meta), path("*.histo"), emit: histo
        path "versions.yml"           ,   emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        jellyfish count \\
            $args \\
            -C \\
            -t $task.cpus \\
            <(zcat ${reads[0]}) <(zcat ${reads[1]}) \\
            -o ${prefix}.jf

        jellyfish histo \\
            -t $task.cpus \\
            ptilo.jf > ${prefix}.histo

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            jellyfish: \$(jellyfish --version 2>&1 | sed 's/^.*jellyfish //' )
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.jf
        touch ${prefix}.histo

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            jellyfish: \$(jellyfish --version 2>&1 | sed 's/^.*jellyfish //' )
        END_VERSIONS
    """
}
