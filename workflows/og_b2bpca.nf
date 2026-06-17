/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { B2B_PREDICTIONS     } from '../subworkflows/local/b2b_predictions/main'
include { ALIGNMENT_DISTANCES } from '../subworkflows/local/alignment_distances/main'
include { PCA_ANALYSIS        } from '../subworkflows/local/pca_analysis/main'
include { SEQUENCE_IDENTITY_ANALYSIS } from '../subworkflows/local/sequence_identity_analysis/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OG_B2BPCA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OG_B2BPCA {

    main:
    ch_versions = Channel.empty()

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
    // Filter OGs with fewer than params.min_seqs sequences
    //
    ch_filtered = ch_og_fastas
        .filter { meta, fasta ->
            def count = fasta.countFasta()
            if (count < params.min_seqs) {
                log.warn "Skipping ${meta.id}: ${count} sequences (< ${params.min_seqs})"
                return false
            }
            return true
        }

    //
    // SUBWORKFLOW: compute b2bTools predictions across all OGs in global batches
    //

    B2B_PREDICTIONS(ch_filtered, params.og_dir)
    ch_versions = ch_versions.mix(B2B_PREDICTIONS.out.versions)

    //
    // SUBWORKFLOW: MAFFT global MSA + MMseqs2 pairwise alignment → dual b2b distance matrices
    //
    ALIGNMENT_DISTANCES(
        ch_filtered,
        B2B_PREDICTIONS.out.og_b2b
    )
    ch_versions = ch_versions.mix(ALIGNMENT_DISTANCES.out.versions)

   //
    // SUBWORKFLOW: PCA + silhouette clustering + summary bar chart
    //
    PCA_ANALYSIS(ALIGNMENT_DISTANCES.out.distance_matrices)
    ch_versions = ch_versions.mix(PCA_ANALYSIS.out.versions)
    //
    // SUBWORKFLOW: Sequence identity matrices and heatmaps from MAFFT and MMseqs2
    // This is a separate subworkflow that reuses MAFFT and MMseqs2 outputs
    // It is b2b backbone independent and provides additional information for data understanding
    // Note: Uses distance matrix sequence ordering and cluster labels to organize heatmaps consistently
    //
    SEQUENCE_IDENTITY_ANALYSIS(
        ALIGNMENT_DISTANCES.out.mafft_alignments,
        ALIGNMENT_DISTANCES.out.mmseqs_pairali,
        PCA_ANALYSIS.out.cluster_labels,
        ALIGNMENT_DISTANCES.out.distance_matrices
    )
    ch_versions = ch_versions.mix(SEQUENCE_IDENTITY_ANALYSIS.out.versions)

    emit:
    versions = ch_versions
    identity_matrices = SEQUENCE_IDENTITY_ANALYSIS.out.identity_matrices
    identity_heatmaps = SEQUENCE_IDENTITY_ANALYSIS.out.identity_heatmaps
}
