/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B2BTOOLS — run b2bTools predictions on a batch FASTA
    ─────────────────────────────────────────────────────────────────────────────────
    Predictors: DynaMine, DisoMine, EFoldMine, AgMata  (PSPer excluded)
    Output columns: sequence_id, residue, residue_index, backbone, sidechain, helix,
                    coil, sheet, ppII, disoMine, earlyFolding, agmata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process B2BTOOLS {
    tag "${meta.id}"
    label 'process_medium'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_b2b.tsv"), emit: predictions
    path 'versions.yml'                         , emit: versions

    script:
    """
    b2bTools \\
        -i ${fasta} \\
        -o ${meta.id}.json \\
        -t ${meta.id}_b2b.tsv \\
        --sep tab \\
        --dynamine \\
        --disomine \\
        --efoldmine \\
        --agmata

    rm -f ${meta.id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        b2bTools: \$(b2bTools --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_b2b.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        b2bTools: 3.0.8
    END_VERSIONS
    """
}
