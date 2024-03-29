process SPLITREADS {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
    //               https://github.com/nf-core/modules/tree/master/modules
    //               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
    //               https://nf-co.re/join
    // TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
    //               All other parameters MUST be provided using the "task.ext" directive, see here:
    //               https://www.nextflow.io/docs/latest/process.html#ext
    //               where "task.ext" is a string.
    //               Any parameters that need to be evaluated in the context of a particular sample
    //               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
    // TODO nf-core: Software that can be piped together SHOULD be added to separate module files
    //               unless there is a run-time, storage advantage in implementing in this way
    //               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
    //                 bwa mem | samtools view -B -T ref.fasta
    // TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
    //               list (`[]`) instead of a file can be used to work around this issue.

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // cv is a vanilla container - nothing installed
    conda "bioconda::seqkit=2.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(reads), path(pairedreads), path(unpairedreads)

    output:
        tuple val(meta), path("output/pairused*.fastq.gz")  , emit: usedreadsp
        tuple val(meta), path("output/unpairused*.fastq.gz"), emit: usedreadsup
        tuple val(meta), path("output/unused*.fastq.gz")    , emit: unusedreads
        path "versions.yml"                                 , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = '0.0.1'
    """
    SEQKIT_THREADS=$task.cpus

    mkdir output

    # Create paths with standardised names
    [ ! -f ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

    [ ! -f output/pairused.${prefix}_1.fastq.gz ] && ln -sf \$(readlink -f ${pairedreads[0]}) output/pairused.${prefix}_1.fastq.gz
    [ ! -f output/pairused.${prefix}_2.fastq.gz ] && ln -sf \$(readlink -f ${pairedreads[1]}) output/pairused.${prefix}_2.fastq.gz

    # Extract headers of paired and unpaired reads
    for f in $pairedreads; do tar -Oxzf \$f | seqkit fx2tab -n; done | sed -e "s/[ \\/].*\$//" | sort | uniq > pairedhead.txt
    for f in $unpairedreads; do tar -Oxzf \$f | seqkit fx2tab -n; done | sed -e "s/[ \\/].*\$//" | sort | uniq > unpairedhead.txt
    cat pairedhead.txt unpairedhead.txt > usedhead.txt

    # Extract reads
    for i in 1 2; do
        zcat ${prefix}_\${i}.fastq.gz | seqkit grep -f unpairedhead.txt | gzip > "unpairused.${prefix}_\${i}.fastq.gz"
        zcat ${prefix}_\${i}.fastq.gz | seqkit grep -v -f usedhead.txt | gzip > "unused.${prefix}_\${i}.fastq.gz"
    done

    # Check pairing

    regex=\$(if [ \$( zcat ${reads[0]} | head -1 | grep -e "\\/[12]") ]; then echo "--id-regexp '^(\\S+)\\/[12]'"; else echo ""; fi)
    for t in unpairused unused
    do
        if [[ -s "\$t.${prefix}_1.fastq.gz" || -s "\$t.${prefix}_2.fastq.gz" ]]
        then
            for i in 1 2; do cp \$t.${prefix}_\${i}.fastq.gz output/\$t.${prefix}_\${i}.fastq.gz; done
        else
            seqkit pair \\
                -1 \$t.${prefix}_1.fastq.gz \\
                -2 \$t.${prefix}_2.fastq.gz \\
                -O output/ \\
                \$regex \\
                --force \\
                2> \$t.pairing.log
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit version 2>&1) )
    END_VERSIONS

    """

}
