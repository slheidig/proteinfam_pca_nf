/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEQUENCE_IDENTITY_MATRIX — compute pairwise sequence identity matrix and heatmap
    ─────────────────────────────────────────────────────────────────────────────────
    For MAFFT: Calculate sequence identity from MSA alignment
    For MMseqs2: Read pairwise sequence identity from pident column
    
    Output:
      - {og_id}_mafft_identity_matrix.csv / {og_id}_mafft_identity_heatmap.png
      - {og_id}_mmseqs2_identity_matrix.csv / {og_id}_mmseqs2_identity_heatmap.png
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SEQUENCE_IDENTITY_MATRIX {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(ali_file), val(mode), path(sequence_identity_matrix_script)
    path cluster_csv, optional: true
    path ref_matrix, optional: true
    path sequence_order, optional: true

    output:
    tuple val(meta), path("${meta.id}_${mode}_identity_matrix.csv")  , emit: matrix
    tuple val(meta), path("${meta.id}_${mode}_identity_heatmap.png"), optional: true, emit: heatmap
    path 'versions.yml'                                               , emit: versions

    script:
    def cluster_arg = cluster_csv ? "--cluster-labels ${cluster_csv}" : ""
    def ref_arg = ref_matrix ? "--ref-matrix ${ref_matrix}" : ""
    def seq_order_arg = sequence_order ? "--sequence-order ${sequence_order}" : ""
    
    """
    
    python3 ${sequence_identity_matrix_script} \\
        --ali         ${ali_file} \\
        --mode        ${mode} \\
        --og-id       ${meta.id} \\
        --out-matrix  ${meta.id}_${mode}_identity_matrix.csv \\
        --out-heatmap ${meta.id}_${mode}_identity_heatmap.png \
        ${cluster_arg} \
        ${ref_arg} \
        ${seq_order_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${mode}_identity_matrix.csv
    touch ${meta.id}_${mode}_identity_heatmap.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        numpy: 1.26.4
        pandas: 2.2.3
        matplotlib: 3.9.4
    END_VERSIONS
    """
}
