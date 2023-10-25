process CATFASTA {
    label 'process_low'

    conda "conda-forge::gzip=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        path(files)

    output:
        path("concat.fasta.gz"), emit: path

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
        """
        cat $files > concat.fasta.gz
        """

}
