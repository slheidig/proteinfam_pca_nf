/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { B2B_PREDICTIONS     } from '../subworkflows/local/b2b_predictions/main'
include { ALIGNMENT_DISTANCES } from '../subworkflows/local/alignment_distances/main'
include { PCA_ANALYSIS        } from '../subworkflows/local/pca_analysis/main'

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
        // DEBUG: Inspect what passes the filter
        .map { meta, fasta ->
            log.info "[FILTERED] OG: ${meta.id} | File: ${fasta} | Exists: ${fasta.exists()}"
            [meta, fasta]
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
        B2B_PREDICTIONS.out.og_b2b,
        params.chunk_size
    )
    ch_versions = ch_versions.mix(ALIGNMENT_DISTANCES.out.versions)

    //
    // SUBWORKFLOW: PCA + silhouette clustering + summary bar chart
    //
    PCA_ANALYSIS(ALIGNMENT_DISTANCES.out.distance_matrices)
    ch_versions = ch_versions.mix(PCA_ANALYSIS.out.versions)

    emit:
    versions = ch_versions
}
