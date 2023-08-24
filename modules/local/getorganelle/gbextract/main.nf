process GBEXTRACT {
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
    conda "bioconda::getorganelle=1.7.7.0 conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        path infile
        //val version //TODO add in a variable specifying the version number of databases to get

    output:
        path("sequences.fasta.gz"),           emit: seqs
        path("genes.cds/gene/gene.fasta.gz"), emit: genes

    when:
        task.ext.when == null || task.ext.when

    script:
    // https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference
    def args = task.ext.args ?: ''
        """
        # Convert the genbank to a fasta
        cmd="import sys; import Bio.SeqIO as io; io.convert('$infile', 'genbank', sys.stdout, 'fasta')"
        python -c "\${cmd}" | gzip > sequences.fasta.gz

        # Extract the genes
        get_annotated_regions_from_gb.py $infile -o genes.cds -t CDS --mix

        gzip -n genes.cds/gene/gene.fasta

            cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            getorganelle: \$(get_organelle_from_reads.py -v 2>&1 | sed -e "s/^.* v//g")
        END_VERSIONS
        """
}
