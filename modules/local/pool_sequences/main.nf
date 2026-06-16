/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POOL_SEQUENCES — tag OG headers and concatenate all FASTAs into one file
    ─────────────────────────────────────────────────────────────────────────────────
    Input  : path to directory containing OG FASTA files (*.fa)
    Output : pooled.fa  — single FASTA; headers have format  >OG_ID|original_seq_id
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process POOL_SEQUENCES {
    label 'process_medium'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path og_dir

    output:
    path 'pooled.fa'    , emit: pooled_fasta
    path 'pooled.names.tsv', emit: names_map
    path 'versions.yml' , emit: versions

    script:
    // pool_sequences.py receives all *.fa files from the input directory
    """
    pool_sequences.py ${og_dir}/*.fa pooled.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch pooled.fa
    touch pooled.names.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
    END_VERSIONS
    """
}
