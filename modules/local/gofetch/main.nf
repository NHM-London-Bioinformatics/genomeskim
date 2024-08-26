process GOFETCH {
    tag "${taxon}"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.1--pyhdfd78af_0':
        'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0' }"

    input:
        val taxon

    output:
        path "output/database/gene.fasta.gz", emit: genes
        path "output/database/seed.fasta.gz", emit: seqs
        // TODO add version if go_fetch has one
        //path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """

        fullargs="--download --output output/ --db database $args"

        python go_fetch.py \\
        --taxonomy $taxon \\
        \$fullargs

        gzip output/database/*.fasta
        gzip output/database/*.gb

        #cat <<-END_VERSIONS > versions.yml
        #"${task.process}":
        #    gofetch: \$(echo \$(go_fetch.py --version 2>&1) ))
        #END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        """

        mkdir -p output/database/

        touch output/database/gene.fasta.gz
        touch output/database/seed.fasta.gz

        #cat <<-END_VERSIONS > versions.yml
        #"${task.process}":
        #    gofetch: \$(echo \$(go_fetch.py --version 2>&1) ))
        #END_VERSIONS
        """

}
