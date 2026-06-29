/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUMMARY_PLOT — aggregate per-OG cluster outputs into bar charts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SUMMARY_PLOT {
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path cluster_csvs

    output:
    path 'cluster_counts_summary.csv', emit: summary_csv
    path 'cluster_count_histogram.csv', emit: histogram_csv
    path 'cluster_counts_barplot.png', optional: true, emit: summary_plot
    path 'versions.yml'              , emit: versions

    script:
    """
    summary_plot.py \\
        --clusters ${cluster_csvs} \\
        --out-csv  cluster_counts_summary.csv \\
        --out-hist-csv cluster_count_histogram.csv \\
        --out-plot cluster_counts_barplot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch cluster_counts_summary.csv
    touch cluster_count_histogram.csv
    touch cluster_counts_barplot.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        pandas: 2.2.3
        matplotlib: 3.9.4
    END_VERSIONS
    """
}
