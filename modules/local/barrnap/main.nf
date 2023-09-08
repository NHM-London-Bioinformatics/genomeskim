
process BARRNAP {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::barrnap=0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/barrnap:0.9--hdfd78af_4' :
        'biocontainers/barrnap:0.9--hdfd78af_4' }"

    input:
        val(meta), path(contigs)

    output:
        path( "*.matches.txt" ) , emit: matches
        path( "*rrna.gff" )     , emit: gff
        path "versions.yml"     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        barrnap \\
            --threads $task.cpus \\
            --outseq ${prefix}_barrnap.fasta \\
            $args \\
            < $contigs \\
            1> ${prefix}_barrnap_rrna.gff \\
            2> barrnap.log

        #this fails when processing an empty file, so it requires a workaround!
        #based on barrnap from ampliseq
        if [ -s ${prefix}_barrnap.fasta ]; then
            grep -h '>' ${prefix}_barrnap.fasta | sed 's/^>//' | sed 's/:\\+/\\t/g' | awk '{print \$2}' | sort -u >\${KINGDOM}.matches.txt
        else
            touch ${prefix}_barrnap_matches.txt
        fi


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            barrnap: \$(echo \$(barrnap --version 2>&1) | sed "s/^.*barrnap //g")
        END_VERSIONS
    """
}
