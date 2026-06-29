#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    og-b2bpca
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Detects structural-dynamic clusters within orthogroups using b2bTools
    biophysical predictions, global MSA (MAFFT), and pairwise alignment (MMseqs2).
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.outdir = "results"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW - All processes consolidated into one file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// =============================================================================
// PROCESSES
// =============================================================================

// POOL_SEQUENCES — tag OG headers and concatenate all FASTAs into one file
process POOL_SEQUENCES {
    label 'process_single'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path og_dir

    output:
    path 'pooled.fa'    , emit: pooled_fasta
    path 'pooled.names.tsv', emit: names_map

    script:
    // pool_sequences.py receives all *.fa files from the input directory
    """
    pool_sequences.py ${og_dir}/*.fa pooled.fa
    """

    stub:
    """
    touch pooled.fa
    touch pooled.names.tsv
    """
}

// B2BTOOLS — run b2bTools predictions on a batch FASTA
process B2BTOOLS {
    tag "${meta.id}"
    label 'process_single'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_b2b.tsv"), emit: predictions

    script:
    """
    b2bTools \\
        -i ${fasta} \\
        -o ${meta.id}.json \\
        -t ${meta.id}_b2b.tsv \\
        --sep tab \\
        --dynamine \\
        --disomine \\
        --efoldmine
    rm -f ${meta.id}.json
    """

    stub:
    """
    touch ${meta.id}_b2b.tsv
    """
}

// REATTACH_B2B — split the merged b2bTools TSV back into per-OG files
process REATTACH_B2B {
    label 'process_highmem'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path merged_tsv
    path names_map

    output:
    path '*_b2b.tsv', emit: og_predictions

    script:
    """
    reattach_b2b.py ${merged_tsv} ${names_map}
    """

    stub:
    """
    touch stub_b2b.tsv
    """
}

// MAFFT — global multiple sequence alignment with --reorder
process MAFFT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mafft:7.525--h031d066_1'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln.fa"), emit: alignment

    script:
    def args = task.ext.args ?: '--reorder --auto'
    def og_id = fasta.baseName
    """
    mafft ${args} --thread ${task.cpus} "${fasta}" > "${og_id}.aln.fa"
    """

    stub:
    """
    touch ${meta.id}.aln.fa
    """
}

// MMSEQS2_EASYSEARCH — all-vs-all pairwise alignment within an OG
process MMSEQS2_EASYSEARCH {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.pairali.tsv"), emit: pairali

    script:
    def args = task.ext.args ?: ''
    def og_id = fasta.baseName
    """
    mkdir -p tmp
    mmseqs easy-search \\
        "${fasta}" \\
        "${fasta}" \\
        "${og_id}.pairali.tsv" \\
        tmp \\
        --format-output "query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,pident,evalue" \\
        --threads ${task.cpus} \\
        ${args}
    rm -rf tmp
    """

    stub:
    """
    touch ${meta.id}.pairali.tsv
    """
}

// B2B_DISTANCE_MATRIX — build NxN b2b backbone distance matrix for one OG
process B2B_DISTANCE_MATRIX {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(b2b_tsv), path(ali_file), val(mode)

    output:
    tuple val(meta), path("${meta.id}_${mode}_b2b_dist.csv"), emit: matrix
    tuple val(meta), path("${meta.id}_${mode}_pair_meta.tsv"), emit: pair_meta

    script:
    """
    b2b_distance_matrix.py \\
        --b2b        ${b2b_tsv} \\
        --ali        ${ali_file} \\
        --mode       ${mode} \\
        --og-id      ${meta.id} \\
        --out-matrix ${meta.id}_${mode}_b2b_dist.csv \\
        --out-meta   ${meta.id}_${mode}_pair_meta.tsv
    """

    stub:
    """
    touch ${meta.id}_${mode}_b2b_dist.csv
    touch ${meta.id}_${mode}_pair_meta.tsv
    """
}

// SEQUENCE_IDENTITY_MATRIX — compute pairwise sequence identity matrix only
process SEQUENCE_IDENTITY_MATRIX {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(ali_file), val(mode), path(sequence_identity_matrix_script)

    output:
    tuple val(meta), path("${meta.id}_${mode}_identity_matrix.csv"), emit: matrix

    script:
    """
    python3 ${sequence_identity_matrix_script} \\
        --ali         ${ali_file} \\
        --mode        ${mode} \\
        --og-id       ${meta.id} \\
        --out-matrix  ${meta.id}_${mode}_identity_matrix.csv
    """

    stub:
    """
    touch ${meta.id}_${mode}_identity_matrix.csv
    """
}

// PCA_CLUSTERING — PCA + silhouette-based KMeans clustering per matrix
process PCA_CLUSTERING {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(matrix_csv)

    output:
    tuple val(meta), path("plots/${meta.id}_${meta.mode}_pca.png"), optional: true, emit: pca_plot
    tuple val(meta), path("csv/${meta.id}_${meta.mode}_clusters.csv"), emit: cluster_labels
    tuple val(meta), path("csv/${meta.id}_${meta.mode}_pca_meta.csv"), emit: pca_meta
    tuple val(meta), path("csv/${meta.id}_${meta.mode}_sequence_order.txt"), emit: sequence_order

    def external_labels_arg = params.external_labels ? "--external-labels ${params.external_labels}" : ""

    script:
    """
    mkdir -p plots csv
    pca_clustering.py \\
        --matrix       ${matrix_csv} \\
        --og-id        ${meta.id} \\
        --mode         ${meta.mode} \\
        --out-plot     "plots/${meta.id}_${meta.mode}_pca.png" \\
        --out-clusters "csv/${meta.id}_${meta.mode}_clusters.csv" \\
        --out-meta     "csv/${meta.id}_${meta.mode}_pca_meta.csv" \\
        --out-sequence-order "csv/${meta.id}_${meta.mode}_sequence_order.txt" \
        ${external_labels_arg} \
    """

    stub:
    """
    mkdir -p plots csv
    touch plots/${meta.id}_${meta.mode}_pca.png
    touch csv/${meta.id}_${meta.mode}_clusters.csv
    touch csv/${meta.id}_${meta.mode}_pca_meta.csv
    touch csv/${meta.id}_${meta.mode}_sequence_order.txt
    """
}

// PLOT_BOTH_HEATMAPS — plot both b2b distance and sequence identity heatmaps
process PLOT_BOTH_HEATMAPS {
    tag "${meta.id} [${meta.mode}]"
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    tuple val(meta), path(b2b_matrix), path(seq_matrix), path(cluster_csv), path(seq_order)

    output:
    tuple val(meta), path("plots/${meta.id}_${meta.mode}_b2b_distance_heatmap.png"), emit: b2b_heatmap
    tuple val(meta), path("plots/${meta.id}_${meta.mode}_seq_identity_heatmap.png"), emit: seq_heatmap

    script:
    """
    mkdir -p plots
    plot_both_heatmaps.py \
        --b2b-matrix ${b2b_matrix} \
        --seq-id-matrix ${seq_matrix} \
        --cluster-labels ${cluster_csv} \
        --sequence-order ${seq_order} \
        --og-id ${meta.id} \
        --mode ${meta.mode} \
        --out-b2b-heatmap "plots/${meta.id}_${meta.mode}_b2b_distance_heatmap.png" \
        --out-seq-heatmap "plots/${meta.id}_${meta.mode}_seq_identity_heatmap.png"
    """

    stub:
    """
    mkdir -p plots
    touch plots/${meta.id}_${meta.mode}_b2b_distance_heatmap.png
    touch plots/${meta.id}_${meta.mode}_seq_identity_heatmap.png
    """
}

// SUMMARY_PLOT — aggregate per-OG cluster outputs into bar charts
process SUMMARY_PLOT {
    label 'process_low'
    container 'docker.io/slheidig/og_b2b_pca:latest'

    input:
    path cluster_csvs

    output:
    path "cluster_counts_summary.csv", emit: summary_csv
    path "cluster_count_histogram.csv", emit: histogram_csv
    path "cluster_counts_barplot.png", optional: true, emit: summary_plot

    script:
    """
    summary_plot.py \\
        --clusters ${cluster_csvs} \\
        --out-csv  cluster_counts_summary.csv \\
        --out-hist-csv cluster_count_histogram.csv \\
        --out-plot cluster_counts_barplot.png
    """

    stub:
    """
    touch cluster_counts_summary.csv
    touch cluster_count_histogram.csv
    touch cluster_counts_barplot.png
    """
}

// =============================================================================
// MAIN WORKFLOW
// =============================================================================

workflow OG_B2BPCA {

    main:
    def og_dir = file(params.og_dir)
    def fa_files = og_dir.list().findAll { it.toString().endsWith('.fa') }.collect { og_dir / it }
    log.info "Creating channel from ${fa_files.size()} .fa files"

    ch_og_fastas = Channel.from(fa_files)
        .map { fasta ->
            def og_id = fasta.baseName
            def meta  = [id: og_id]
            [meta, fasta]
        }
   
    //
    // Filter OGs with fewer than params.min_seqs sequences (skip if aligned)
    //
    if (params.aligned) {
        ch_filtered = ch_og_fastas
        log.info "Skipping sequence filtering (aligned=true)"
    } else {
        ch_filtered = ch_og_fastas
            .filter { meta, fasta ->
                def count = fasta.countFasta()
                if (count < params.min_seqs) {
                    log.warn "Skipping ${meta.id}: ${count} sequences (< ${params.min_seqs})"
                    return false
                }
                return true
            }
    }

    //
    // STAGE 1: B2B PREDICTIONS
    // Pool all filtered OG FASTAs, split into batches, run b2bTools, reattach to OGs
    //
    
    // Pool sequences
    POOL_SEQUENCES(og_dir)
    ch_names = POOL_SEQUENCES.out.names_map

    // Split into batches
    ch_batches = POOL_SEQUENCES.out.pooled_fasta
        .splitFasta(by: params.b2b_batch_size, file: true)
        .map { batch_file ->
            def batch_id = "batch_${batch_file.baseName}"
            [ [id: batch_id], batch_file ]
        }

    // Run b2bTools on batches
    B2BTOOLS(ch_batches)

    // Merge and reattach
    ch_merged = B2BTOOLS.out.predictions
        .map { meta, tsv -> tsv }
        .collectFile(name: 'all_b2b_predictions.tsv', keepHeader: true, skip: 1)

    REATTACH_B2B(ch_merged, ch_names)

    // Reconstruct per-OG files
    ch_og_b2b = REATTACH_B2B.out.og_predictions
        .flatten()
        .map { tsv_file ->
            def og_id = tsv_file.baseName.replace('_b2b', '')
            [ [id: og_id], tsv_file ]
        }

    //
    // STAGE 2: ALIGNMENT DISTANCES
    // MAFFT global MSA + MMseqs2 pairwise alignment
    //
    def script_file = file("${workflow.projectDir}/bin/sequence_identity_matrix.py")

    if (params.aligned) {
        ch_mafft_out = ch_og_fastas.map { meta, fasta -> [meta, fasta] }
        ch_mmseqs_out = ch_og_fastas.map { meta, fasta -> [meta, fasta] }
    } else {
        // MAFFT alignment
        MAFFT(ch_filtered)
        ch_mafft_out = MAFFT.out.alignment
            .map { meta, aln_file -> [ meta, aln_file ] }

        // MMseqs2 pairwise alignment
        MMSEQS2_EASYSEARCH(ch_filtered)
        ch_mmseqs_out = MMSEQS2_EASYSEARCH.out.pairali
            .map { meta, pairali_file -> [ meta, pairali_file ] }
    }

    // Join b2b with alignments and build distance matrices
    ch_b2b_mafft = ch_og_b2b
        .join(ch_mafft_out)
        .map { meta, b2b_tsv, aln_fa ->
            [ meta.plus([mode: 'mafft']), b2b_tsv, aln_fa, 'mafft' ]
        }

    ch_b2b_mmseqs = ch_og_b2b
        .join(ch_mmseqs_out)
        .map { meta, b2b_tsv, aln_fa ->
            [ meta.plus([mode: 'mmseqs2']), b2b_tsv, aln_fa, 'mmseqs2' ]
        }

    B2B_DISTANCE_MATRIX(ch_b2b_mafft.mix(ch_b2b_mmseqs))

    // Compute sequence identity matrices
    ch_mafft_with_mode = ch_mafft_out.map { meta, aln_file -> [ meta.plus([mode: 'mafft']), aln_file ] }
    ch_mmseqs_with_mode = ch_mmseqs_out.map { meta, pairali_file -> [ meta.plus([mode: 'mmseqs2']), pairali_file ] }
    
    ch_identity_inputs = ch_mafft_with_mode
        .map { meta, aln_file -> [ meta, aln_file, meta.mode, script_file ] }
        .mix(
            ch_mmseqs_with_mode
                .map { meta, pairali_file -> [ meta, pairali_file, meta.mode, script_file ] }
        )
    
    SEQUENCE_IDENTITY_MATRIX(ch_identity_inputs)
    
    // Store identity matrices for later heatmap generation
    ch_identity_matrices_all = SEQUENCE_IDENTITY_MATRIX.out.matrix

    //
    // STAGE 3: PCA ANALYSIS
    // PCA + clustering on distance matrices
    //
    if (params.distance_type == 'seq_distance') {
        ch_matrices_for_pca = SEQUENCE_IDENTITY_MATRIX.out.matrix
    } else {
        ch_matrices_for_pca = B2B_DISTANCE_MATRIX.out.matrix
    }

    ch_matrices_for_pca | PCA_CLUSTERING

    // Collect cluster labels for summary plot
    SUMMARY_PLOT(
        PCA_CLUSTERING.out.cluster_labels
            .map { meta, csv -> csv }
            .collect()
    )
    
    // Generate both heatmaps with cluster ordering
    // Join b2b matrices, seq identity matrices, cluster labels, and sequence order
    ch_heatmap_inputs = B2B_DISTANCE_MATRIX.out.matrix
        .join(SEQUENCE_IDENTITY_MATRIX.out.matrix)
        .join(PCA_CLUSTERING.out.cluster_labels)
        .join(PCA_CLUSTERING.out.sequence_order)
        .map { meta, b2b_csv, seq_csv, cluster_csv, seq_order ->
            [ meta, b2b_csv, seq_csv, cluster_csv, seq_order ]
        }
    
    PLOT_BOTH_HEATMAPS(ch_heatmap_inputs)
    ch_b2b_heatmaps = PLOT_BOTH_HEATMAPS.out.b2b_heatmap
    ch_seq_heatmaps = PLOT_BOTH_HEATMAPS.out.seq_heatmap
    
    // Combine heatmaps for downstream use
    ch_identity_heatmaps = ch_seq_heatmaps

    //
    // STAGE 4: SEQUENCE IDENTITY ANALYSIS
    //
    if (params.distance_type == 'seq_distance') {
        SEQUENCE_IDENTITY_ANALYSIS_EXPANDED(
            ch_mafft_out,
            ch_mmseqs_out,
            PCA_CLUSTERING.out.cluster_labels,
            B2B_DISTANCE_MATRIX.out.matrix,
            PCA_CLUSTERING.out.sequence_order,
            ch_identity_matrices_all,
            ch_identity_heatmaps
        )
    } else {
        SEQUENCE_IDENTITY_ANALYSIS_EXPANDED(
            ch_mafft_out,
            ch_mmseqs_out,
            PCA_CLUSTERING.out.cluster_labels,
            B2B_DISTANCE_MATRIX.out.matrix,
            PCA_CLUSTERING.out.sequence_order,
            Channel.empty(),
            ch_identity_heatmaps
        )
    }

    emit:
    identity_matrices = SEQUENCE_IDENTITY_ANALYSIS_EXPANDED.out.identity_matrices
    identity_heatmaps = SEQUENCE_IDENTITY_ANALYSIS_EXPANDED.out.identity_heatmaps
}

// =============================================================================
// SEQUENCE_IDENTITY_ANALYSIS_EXPANDED - Inlined version of the subworkflow
// =============================================================================

workflow SEQUENCE_IDENTITY_ANALYSIS_EXPANDED {

    def script_file = file("${workflow.projectDir}/bin/sequence_identity_matrix.py")

    take:
    ch_mafft_alignments
    ch_mmseqs_pairali
    ch_cluster_labels
    ch_distance_matrices
    ch_sequence_order
    ch_precomputed_matrices
    ch_precomputed_heatmaps

    main:
    // Use pre-computed identity matrices and heatmaps
    ch_identity_matrices = ch_precomputed_matrices
    ch_identity_heatmaps = ch_precomputed_heatmaps

    emit:
    identity_matrices = ch_identity_matrices
    identity_heatmaps  = ch_identity_heatmaps
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    OG_B2BPCA()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
