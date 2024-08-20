process CHECKFORMAT {
    tag "$file"
    label 'process_low'

    conda "conda-forge::gzip=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        path(file)

    output:
        env type         , emit: format
        path('type.txt') , emit: formatfile

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        type=""

        filepath=$file
        if [[ "$file" == *".gz" ]]
        then
            gzip -c -d "$file" > "\${filepath%.gz}"
            filepath="\${filepath%.gz}"
        fi

        if \$(head -1 "\$filepath" | grep -q "^>")
        then
            type="fasta"
        elif \$(head -1 "\$filepath" | grep -q "^>")
        then
            type="fastq"
        elif \$(head -1 "\$filepath" | grep -q "^LOCUS")
        then
            type="gb"
        else
            type="unknown"
        fi

        if [[ "$file" == *".gz" ]]
        then
            rm "\$filepath"
        fi

        echo "\$type" > type.txt

        """

    stub:
        def args = task.ext.args ?: ''
        """
        type="unknown"

        echo "\$type" > type.txt

        """
}
