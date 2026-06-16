/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REATTACH_B2B — split the merged b2bTools TSV back into per-OG files
    ─────────────────────────────────────────────────────────────────────────────────
    Input : merged TSV of all batch predictions (sequence_id = OG_ID|original_seq_id)
    Output: one *_b2b.tsv per OG; sequence_id column restored to original IDs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REATTACH_B2B {
    label 'process_medium'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path merged_tsv
    path names_map

    output:
    path '*_b2b.tsv'    , emit: og_predictions
    path 'versions.yml' , emit: versions

    script:
    """
    reattach_b2b.py ${merged_tsv} ${names_map}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch stub_b2b.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
    END_VERSIONS
    """
}
