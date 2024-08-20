process CATREADS {
    tag "$meta.id"
    label 'process_low'

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

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
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

        grep "_1" readfiles.txt > r1.txt
        grep "_2" readfiles.txt > r2.txt

        if [ \$(cat r1.txt | wc -l) != \$(cat r2.txt | wc -l) ]
        then
            echo "ERROR: not the same number of files for each read direction"
            exit 1
        fi

        touch "output/${prefix}_\${outname}_1.\${ext}.gz"
        touch "output/${prefix}_\${outname}_2.\${ext}.gz"
        """
}
