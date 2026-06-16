/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B2B_DISTANCE_MATRIX — build NxN b2b backbone distance matrix for one OG
    ─────────────────────────────────────────────────────────────────────────────────
    mode = 'mafft'   → ali_file is a multiple sequence alignment FASTA
    mode = 'mmseqs2' → ali_file is the MMseqs2 pairali.tsv
    For each unique pair (a, b) from itertools.combinations:
      - map b2b backbone values onto aligned columns
      - skip gap positions
      - distance = Σ|backbone_a − backbone_b| / n_compared_residues
    Tracks: seq_a_len, seq_b_len, ali_len, n_compared per pair.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process B2B_DISTANCE_MATRIX {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(b2b_tsv), path(ali_file), val(mode)

    output:
    tuple val(meta), path("${meta.id}_${mode}_b2b_dist.csv"), emit: matrix
    tuple val(meta), path("${meta.id}_${mode}_pair_meta.tsv"), emit: pair_meta
    path 'versions.yml'                                      , emit: versions

    script:
    """
    b2b_distance_matrix.py \\
        --b2b        ${b2b_tsv} \\
        --ali        ${ali_file} \\
        --mode       ${mode} \\
        --og-id      ${meta.id} \\
        --out-matrix ${meta.id}_${mode}_b2b_dist.csv \\
        --out-meta   ${meta.id}_${mode}_pair_meta.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${mode}_b2b_dist.csv
    touch ${meta.id}_${mode}_pair_meta.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        numpy: 1.26.4
        pandas: 2.2.3
    END_VERSIONS
    """
}
