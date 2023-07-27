process PREPDATABASES {
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
    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        val genometype
        path seedspathin
        path genespathin
        //val version //TODO add in a variable specifying the version number of databases to get

    output:
        path("*Seed.fasta"),  emit: seeds
        path("*Label.fasta"), emit: labels

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}" Don't need this - $meta not defined
    //def seedspath = seedspathin ? "${seedspathin}" : ''
    //def genespath = genespathin ? "${seedspathin}" : ''
    def version = '0.0.1'

        """
        repo="Kinggerm/GetOrganelleDB/master"

        declare -A paths
        paths['seed']=$seedspathin # Does this work?
        paths['gene']=$genespathin # If not uncomment lines 45 and 46 and remove in from the end

        declare -A dbs
        dbs['seed']="Seed"
        dbs['gene']="Label"

        for id in seed gene
        do
        if [[ \${paths[\$id]} == '' ]]
        then
            url="https://raw.githubusercontent.com/\$repo/$version/\${dbs[\$id]}Database/${genometype}.fasta"
            if ! wget -q --spider \$url
            then
                >&2 echo "Error: could not find a default database in github repo \$repo, are version $version and organelle genome $genometype both correct?"
                exit 1
            fi
            wget -qO "default_\${dbs[\$id]}.fasta" \$url
        else
            cp \${paths[\$id]} custom_\${dbs[\$id]}.fasta
        fi
        done

        """

}
