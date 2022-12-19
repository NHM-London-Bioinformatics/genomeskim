process GETORGANELLE {
    tag "$meta.id"
    label 'process_high'

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

    // TODO figure out a way to generate some multiQC-appropriate files
    // TODO set up for single-end reads

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::getorganelle=1.7.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(reads)
    path seeds
    path genes

    output:
    // Output complete contigs
    tuple val(meta), path("output/*path_sequence.fasta"), emit: contig
    tuple val(meta), path("getorganelle.log")               , emit: log
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    if ( params.enable_conda ) {
        """
        get_organelle_from_reads.py \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o output \\
            -t $task.cpus \\
            --zip-files \\
            --verbose \\
            -s $seeds \\
            --genes $genes \\
            $args \\
            2>&1 > getorganelle.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            getorganelle: \$(get_organelle_from_reads.py -v 2>&1 | sed -e "s/^.* v//g")
        END_VERSIONS

        """
    } else {
        """
        echo "Only conda is currently supported by the GetOrganelle module"
        exit 1
        """
    }
}
//TODO Need to find a way to output the error message properly to nextflow if there is an error.
//TODO Is it possible to correct what is printed to the terminal if there's an issue?
