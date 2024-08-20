process CATFASTA {
    label 'process_low'

    conda "conda-forge::gzip=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        path(files)
        val prefix

    output:
        path("*concat.fasta.gz"), emit: catfasta

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        cat $files > ${prefix}_concat.fasta.gz

        """

    stub:
        def args = task.ext.args ?: ''
        """
        touch ${prefix}_concat.fasta.gz
        """
}
