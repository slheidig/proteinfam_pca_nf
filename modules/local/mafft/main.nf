/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAFFT — global multiple sequence alignment with --reorder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MAFFT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mafft:7.525--h031d066_1'

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.aln.fa"), emit: alignment
    path 'versions.yml'                        , emit: versions

    script:
    def args = task.ext.args ?: '--reorder --auto'
    """
    for fasta in ${fastas}; do
        base=\$(basename "\$fasta")
        og_id="\${base%.*}"
        mafft ${args} --thread ${task.cpus} "\$fasta" > "\${og_id}.aln.fa"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | grep -oP 'v\\d+\\.\\d+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.aln.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: 7.525
    END_VERSIONS
    """
}
