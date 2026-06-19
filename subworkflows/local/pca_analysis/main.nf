/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PCA_ANALYSIS subworkflow
    ─────────────────────────────────────────────────────────────────────────────────
    1. Run PCA + silhouette-based KMeans clustering on each b2b distance matrix.
    2. Aggregate all per-OG cluster counts and produce a summary bar chart.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PCA_CLUSTERING } from '../../../modules/local/pca_clustering/main'
include { SUMMARY_PLOT   } from '../../../modules/local/summary_plot/main'

workflow PCA_ANALYSIS {

    take:
    ch_distance_matrices    // channel: [ val(meta), path(*_b2b_dist.csv) ]

    main:
    ch_versions = Channel.empty()

    //
    // PCA + clustering — one job per (OG × alignment-mode) combination
    //
    ch_distance_matrices | PCA_CLUSTERING
    ch_versions = ch_versions.mix(PCA_CLUSTERING.out.versions)

    //
    // Collect all per-OG cluster label CSVs and build a summary bar chart.
    // Runs once after all PCA jobs complete.
    //
    SUMMARY_PLOT(
        PCA_CLUSTERING.out.cluster_labels
            .map { meta, csv -> csv }
            .collect()
    )
    ch_versions = ch_versions.mix(SUMMARY_PLOT.out.versions)

    emit:
    pca_plots      = PCA_CLUSTERING.out.pca_plot       // channel: [ val(meta), path(*_pca.png) ]
    heatmaps       = PCA_CLUSTERING.out.heatmap        // channel: [ val(meta), path(*_heatmap.png) ]
    cluster_labels = PCA_CLUSTERING.out.cluster_labels  // channel: [ val(meta), path(*.csv) ]
    sequence_order = PCA_CLUSTERING.out.sequence_order // channel: [ val(meta), path(*_sequence_order.txt) ]
    versions       = ch_versions
}
