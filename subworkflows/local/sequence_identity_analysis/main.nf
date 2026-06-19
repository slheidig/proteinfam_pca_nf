/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEQUENCE_IDENTITY_ANALYSIS subworkflow
    ─────────────────────────────────────────────────────────────────────────────────
    For each OG:
      1. Compute pairwise sequence identity matrix from MAFFT MSA
      2. Compute pairwise sequence identity matrix from MMseqs2 pairwise alignments
      3. Generate heatmaps for both matrices, sorted by cluster labels
    
    This is a separate subworkflow that reuses outputs from MAFFT and MMseqs2.
    It is b2b backbone independent and provides additional information for data understanding.
    
    Output files:
      - {og_id}_mafft_identity_matrix.csv
      - {og_id}_mafft_identity_heatmap.png
      - {og_id}_mmseqs2_identity_matrix.csv
      - {og_id}_mmseqs2_identity_heatmap.png
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQUENCE_IDENTITY_MATRIX } from '../../../modules/local/sequence_identity_matrix/main'

workflow SEQUENCE_IDENTITY_ANALYSIS {

    def script_file = file("${workflow.projectDir}/bin/sequence_identity_matrix.py")

    take:
    ch_mafft_alignments     // channel: [ val(meta), path(*.aln.fa) ] - meta has id only
    ch_mmseqs_pairali        // channel: [ val(meta), path(*.pairali.tsv) ] - meta has id only
    ch_cluster_labels        // channel: [ val(meta), path(*_clusters.csv) ] - meta has id and mode, optional
    ch_distance_matrices    // channel: [ val(meta), path(*_b2b_dist.csv) ] - optional, for sequence ordering
    ch_sequence_order       // channel: [ val(meta), path(*_sequence_order.txt) ] - sequence order from PCA

    main:
    ch_versions = Channel.empty()

    //
    // Add mode='mafft' to MAFFT alignment meta
    //
    ch_mafft_with_mode = ch_mafft_alignments
        .map { meta, aln_file ->
            [ meta.plus([mode: 'mafft']), aln_file ]
        }

    //
    // Add mode='mmseqs2' to MMseqs2 pairali meta
    //
    ch_mmseqs_with_mode = ch_mmseqs_pairali
        .map { meta, pairali_file ->
            [ meta.plus([mode: 'mmseqs2']), pairali_file ]
        }

    //
    // Build a map of distance matrices by og_id and mode for reference ordering
    //
    ch_distance_map = ch_distance_matrices
        .map { meta, dist_file ->
            def key = "${meta.id}:${meta.mode}"
            [ key, dist_file ]
        }
    
    //
    // Build a map of sequence order files by og_id and mode
    //
    ch_sequence_order_map = ch_sequence_order
        .map { meta, seq_order_file ->
            def key = "${meta.id}:${meta.mode}"
            [ key, seq_order_file ]
        }
    
    //
    // If cluster labels are provided, join them with the alignment files
    // The cluster labels channel has meta with og_id and mode
    //
    if (ch_cluster_labels) {
        // Prepare MAFFT for join: create channel of [join_key, aln_file, meta]
        // where join_key is a string "og_id:mode"
        ch_mafft_for_join = ch_mafft_with_mode
            .map { meta, aln_file ->
                def join_key = "${meta.id}:${meta.mode}"
                [ join_key, aln_file, meta ]
            }
        
        // Prepare MMseqs2 for join
        ch_mmseqs_for_join = ch_mmseqs_with_mode
            .map { meta, pairali_file ->
                def join_key = "${meta.id}:${meta.mode}"
                [ join_key, pairali_file, meta ]
            }
        
        // Prepare cluster labels: [join_key, csv]
        ch_cluster_for_join = ch_cluster_labels
            .map { meta, csv ->
                def join_key = "${meta.id}:${meta.mode}"
                [ join_key, csv ]
            }
        
        // Join MAFFT with cluster labels using join_key
        ch_mafft_with_clusters = ch_mafft_for_join
            .join(ch_cluster_for_join)
            .map { join_key, aln_file, meta, cluster_csv ->
                [ meta, aln_file, cluster_csv, meta.mode ]
            }
        
        // Join MMseqs2 with cluster labels using join_key
        ch_mmseqs_with_clusters = ch_mmseqs_for_join
            .join(ch_cluster_for_join)
            .map { join_key, pairali_file, meta, cluster_csv ->
                [ meta, pairali_file, cluster_csv, meta.mode ]
            }
        
        // Mix and run identity matrix computation
        ch_identity_inputs = ch_mafft_with_clusters.mix(ch_mmseqs_with_clusters)
    } else {
        // No cluster labels available, use null for cluster_csv
        ch_identity_inputs = ch_mafft_with_mode
            .map { meta, aln_file -> [ meta, aln_file, null, meta.mode ] }
            .mix(
                ch_mmseqs_with_mode
                    .map { meta, pairali_file -> [ meta, pairali_file, null, meta.mode ] }
            )
    }
    
    //
    // Join identity inputs with distance matrices and sequence order for reference ordering
    //
    ch_identity_with_ref = ch_identity_inputs
        .map { meta, aln_file, cluster_csv, mode ->
            def key = "${meta.id}:${mode}"
            [ key, meta, aln_file, cluster_csv, mode ]
        }
        .join(ch_distance_map)
        .map { key, meta, aln_file, cluster_csv, mode, ref_matrix ->
            [ key, meta, aln_file, cluster_csv, mode, ref_matrix ]
        }
        .join(ch_sequence_order_map, remainder: true)
        .map { key, meta, aln_file, cluster_csv, mode, ref_matrix, seq_order_file ->
            [ meta, aln_file, cluster_csv, mode, script_file, ref_matrix, seq_order_file ]
        }

    //
    // Compute identity matrices and heatmaps
    //
    SEQUENCE_IDENTITY_MATRIX(ch_identity_with_ref)
    ch_versions = ch_versions.mix(SEQUENCE_IDENTITY_MATRIX.out.versions)

    emit:
    identity_matrices = SEQUENCE_IDENTITY_MATRIX.out.matrix     // channel: [ val(meta), path(*_identity_matrix.csv) ]
    identity_heatmaps  = SEQUENCE_IDENTITY_MATRIX.out.heatmap    // channel: [ val(meta), path(*_identity_heatmap.png) ]
    versions           = ch_versions
}
