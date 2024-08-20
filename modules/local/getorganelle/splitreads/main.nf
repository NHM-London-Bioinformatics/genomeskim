process SPLITREADS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.8.1-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0' :
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

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
        """
        SEQKIT_THREADS=$task.cpus

        mkdir output

        # Create paths with standardised names
        [ ! -f ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

        [ ! -f output/pairused.${prefix}_1.fastq.gz ] && ln -sf \$(readlink -f ${pairedreads[0]}) output/pairused.${prefix}_1.fastq.gz
        [ ! -f output/pairused.${prefix}_2.fastq.gz ] && ln -sf \$(readlink -f ${pairedreads[1]}) output/pairused.${prefix}_2.fastq.gz

        # Extract headers of paired and unpaired reads
        #seems getorganelle outputs are sometimes tar - need to write something with tar -Oxzf to deal
        for f in $pairedreads; do zcat \$f | seqkit fx2tab -n; done | sed -e "s/[ \\/].*\$//" | sort | uniq > pairedhead.txt
        for f in $unpairedreads; do zcat \$f | seqkit fx2tab -n; done | sed -e "s/[ \\/].*\$//" | sort | uniq > unpairedhead.txt
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

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir output

        for i in 1 2
        do
            for j in un unpair pair
            do
                touch output/${j}.${prefix}_${i}.fastq.gz
            done
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit version 2>&1) )
        END_VERSIONS

        """


}
