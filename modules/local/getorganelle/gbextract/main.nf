process GBEXTRACT {
    tag ${infile}
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.1--pyhdfd78af_0':
        'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0' }"

    input:
        path infile

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

    script:
        // https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference
        def args = task.ext.args ?: ''
        """

        mkdir -p genes.cds/gene
        touch sequences.fasta.gz
        touch genes.cds/gene/gene.fasta.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            getorganelle: \$(get_organelle_from_reads.py -v 2>&1 | sed -e "s/^.* v//g")
        END_VERSIONS
        """
}
