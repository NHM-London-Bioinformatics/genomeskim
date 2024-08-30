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
        tuple val(meta), path("catoutput/${meta.id}.bed"), emit: bed
        tuple val(meta), path("catoutput/${meta.id}.gff"), emit: gff
        tuple val(meta), path("catoutput/${meta.id}.fas"), emit: fas
        tuple val(meta), path("catoutput/${meta.id}.faa"), emit: faa
        path "versions.yml"                              , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        """
        args=\$(grep -q "circular" $contig && echo "$args" || echo "$args --linear")

        mkdir catoutput
        mkdir mitosoutput

        runmitos.py \
            --input $contig \
            --outdir mitosoutput/ \
            --refdir ./ \
            --refseqver $db \
            \$args \
            &> mitos.log

        dirs='/'
        if [ ! -f mitosoutput/result.mitos ]
        then
            dirs=\$(ls mitosoutput | sed "s:\$:/:")
        fi
        for d in dirs
        do
            f="mitosoutput/\${d}result"
            for ext in bed faa fas geneorder gff mitos seq
            do
                cat \$f.\$ext >> catoutput/${prefix}.\$ext
            done
            cp \$f.png catoutput/${prefix}_\${d%/}.png
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runmitos.py: \$(runmitos.py --version 2>&1 | sed 's/^.*runmitos.py //' )
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        """
        mkdir catoutput

        for i in bed gff fas faa;
        do
            touch catoutput/$prefix.\$i
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runmitos.py: \$(runmitos.py --version 2>&1 | sed 's/^.*runmitos.py //' )
        END_VERSIONS
        """
}
