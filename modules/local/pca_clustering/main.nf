/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PCA_CLUSTERING — PCA + silhouette-based KMeans clustering per matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PCA_CLUSTERING {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(matrix_csv)

    output:
    tuple val(meta), path("${meta.id}_${meta.mode}_pca.png"),     emit: pca_plot
    tuple val(meta), path("${meta.id}_${meta.mode}_heatmap.png"),  emit: heatmap
    tuple val(meta), path("${meta.id}_${meta.mode}_clusters.csv"), emit: cluster_labels
    tuple val(meta), path("${meta.id}_${meta.mode}_pca_meta.csv"), emit: pca_meta
    tuple val(meta), path("${meta.id}_${meta.mode}_sequence_order.txt"), emit: sequence_order
    path 'versions.yml', emit: versions

    def external_labels_arg = params.external_labels ? "--external-labels ${params.external_labels}" : ""

    script:
    """
    pca_clustering.py \\
        --matrix       ${matrix_csv} \\
        --og-id        ${meta.id} \\
        --mode         ${meta.mode} \\
        --out-plot     ${meta.id}_${meta.mode}_pca.png \\
        --out-heatmap  ${meta.id}_${meta.mode}_heatmap.png \\
        --out-clusters ${meta.id}_${meta.mode}_clusters.csv \\
        --out-meta     ${meta.id}_${meta.mode}_pca_meta.csv \\
        --out-sequence-order ${meta.id}_${meta.mode}_sequence_order.txt \
        ${external_labels_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        sklearn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.mode}_pca.png
    touch ${meta.id}_${meta.mode}_heatmap.png
    touch ${meta.id}_${meta.mode}_clusters.csv
    touch ${meta.id}_${meta.mode}_pca_meta.csv
    touch ${meta.id}_${meta.mode}_sequence_order.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        sklearn: 1.6.1
        matplotlib: 3.9.4
    END_VERSIONS
    """
}
