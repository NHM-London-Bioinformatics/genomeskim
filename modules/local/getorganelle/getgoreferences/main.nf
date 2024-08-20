process GETGOREFS {
    label 'process_low'

    conda "conda-forge::gzip=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        val genometype
        //val version //TODO add in a variable specifying the version number of databases to get

    output:
        path("*Seed.fasta.gz") , emit: seeds
        path("*Label.fasta.gz"), emit: labels
        path "versions.yml"    , emit: versions


    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def version = '0.0.1'
        """
        repo="Kinggerm/GetOrganelleDB/master"

        declare -A dbs
        dbs['seed']="Seed"
        dbs['gene']="Label"

        for id in seed gene
        do
            url="https://raw.githubusercontent.com/\$repo/$version/\${dbs[\$id]}Database/${genometype}.fasta"
            if ! wget -q --spider \$url
            then
                >&2 echo "Error: could not find a default database in github repo \$repo, are version $version and organelle genome $genometype both correct?"
                exit 1
            fi
            wget -qO "default_\${dbs[\$id]}.fasta" \$url
        done

        gzip -n default_*.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GetOrganelleDB: $version
        END_VERSIONS

        """

    stub:
        def args = task.ext.args ?: ''
        def version = '0.0.1'
        """

        touch default_Label.fasta.gz
        touch default_Seed.fasta.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GetOrganelleDB: $version
        END_VERSIONS

        """
}
