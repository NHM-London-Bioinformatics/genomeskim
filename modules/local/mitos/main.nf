process MITOS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mitos=2.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.3--pyhdfd78af_0':
        'biocontainers/mitos:2.1.3--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(contig)
        tuple val(dbmeta), path(db)

    output:
        tuple val(meta), path("${meta.id}.bed"), emit: bed
        tuple val(meta), path("${meta.id}.gff"), emit: gff
        tuple val(meta), path("${meta.id}.fas"), emit: fas
        tuple val(meta), path("${meta.id}.faa"), emit: faa
        path "versions.yml"                       , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        """
        args=\$(grep -q "circular" $contig && echo "$args" || echo "$args --linear")

        mkdir output

        runmitos.py \
            --input $contig \
            --outdir output/ \
            --refdir ./ \
            --refseqver $db \
            \$args \
            &> mitos.log

        for f in output/result*; do mv $f ${prefix}.${f##*.}; done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runmitos.py: \$(runmitos.py --version 2>&1 | sed 's/^.*runmitos.py //' )
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        """
        mkdir output

        for i in bed gff fas faa;
        do
            touch $prefix.\$i
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runmitos.py: \$(runmitos.py --version 2>&1 | sed 's/^.*runmitos.py //' )
        END_VERSIONS
        """
}
