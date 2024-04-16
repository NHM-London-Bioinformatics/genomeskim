process MITOS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mitos=2.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.3--pyhdfd78af_0':
        'biocontainers/mitos:2.1.3--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(contig)
        path(db)

    output:
        path "output/result.bed"      , emit: bed
        path "output/result.gff"      , emit: gff
        path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        args=\$(grep -q "circular" $contig && echo "$args" || echo "$args --linear")

        mkdir output

        runmitos.py \
            --input <( zcat $contig )\
            --outdir output/ \
            --refdir ./ \
            --refseqver $db \
            --zip-files
            \$args \
            &> mitos.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runmitos.py: \$(echo \$(runmitos.py --version 2>&1) | sed 's/^.*runmitos.py //' ))
        END_VERSIONS
        """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
