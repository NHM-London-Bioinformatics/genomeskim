process CATREADS {
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
    conda "conda-forge::gzip=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(reads)
        val(outname)

    output:
        tuple val(meta), path("output/*.gz"), emit: catreads

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = '0.0.1'
        """
        mkdir output

        outname=\$( if [ -z "$outname" ]; then echo "catreads"; else echo "$outname"; fi )

        for r in $reads; do echo \$r; done | sort > readfiles.txt

        ext=""
        if [ "\$(grep -c "fastq\\|fq" readfiles.txt)" == "\$(cat readfiles.txt | wc -l)" ]
        then
            ext="fastq"
        elif [ "\$(grep -c "fasta\\|fa" readfiles.txt)" == "\$(cat readfiles.txt | wc -l)" ]
        then
            ext="fasta"
        else
            echo "ERROR: supplied files do not appear to be of a consistent format"
            exit 1
        fi

        #cat readfiles.txt | tee >(awk '/_1/' > r1.txt) >(awk '/_2/' > r2.txt)
        grep "_1" readfiles.txt > r1.txt
        grep "_2" readfiles.txt > r2.txt

        if [ \$(cat r1.txt | wc -l) != \$(cat r2.txt | wc -l) ]
        then
            echo "ERROR: not the same number of files for each read direction"
            exit 1
        fi

        cat r1.txt | xargs cat > "output/${prefix}_\${outname}_1.\${ext}.gz"
        cat r2.txt | xargs cat > "output/${prefix}_\${outname}_2.\${ext}.gz"
        """

}
