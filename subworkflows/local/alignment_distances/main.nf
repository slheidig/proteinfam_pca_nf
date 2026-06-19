/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGNMENT_DISTANCES subworkflow
    ─────────────────────────────────────────────────────────────────────────────────
    For each OG:
      1. Build a global multiple sequence alignment with MAFFT (--reorder).
      2. Run MMseqs2 easy-search all-vs-all for pairwise alignments.
      3. Build a b2b backbone distance matrix using the MAFFT MSA.
      4. Build a b2b backbone distance matrix using the MMseqs2 pairwise alignments.
    Both matrices are emitted for downstream PCA comparison.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAFFT               } from '../../../modules/local/mafft/main'
include { MMSEQS2_EASYSEARCH  } from '../../../modules/local/mmseqs2_easysearch/main'
include { B2B_DISTANCE_MATRIX } from '../../../modules/local/b2b_distance_matrix/main'

workflow ALIGNMENT_DISTANCES {

    take:
    ch_og_fastas    // channel: [ val(meta), path(fasta) ]
    ch_og_b2b       // channel: [ val(meta), path(*_b2b.tsv) ] 

    main:
    ch_versions = Channel.empty()

   //
    // Global MSA — MAFFT with --reorder (run per OG)
    //
    MAFFT(ch_og_fastas)
    ch_versions = ch_versions.mix(MAFFT.out.versions)

    //
    // MAFFT output: preserve original meta for join compatibility
    //
    ch_mafft_out = MAFFT.out.alignment
        .map { meta, aln_file ->
            [ meta, aln_file ]
        }

    //
    // Pairwise alignment — MMseqs2 easy-search all-vs-all per OG
    //
    MMSEQS2_EASYSEARCH(ch_og_fastas)
    ch_versions = ch_versions.mix(MMSEQS2_EASYSEARCH.out.versions)

    //
    // MMSEQS2 output: preserve original meta for join compatibility
    //
    ch_mmseqs_out = MMSEQS2_EASYSEARCH.out.pairali
        .map { meta, pairali_file ->
            [ meta, pairali_file ]
        }


    //
    // Join b2b TSV with MAFFT alignment (key: meta.id)
    // Augment meta with mode so output filenames and PCA labels are unambiguous.
    //
    ch_b2b_mafft = ch_og_b2b
        .join(ch_mafft_out)
        .map { meta, b2b_tsv, aln_fa ->
            [ meta.plus([mode: 'mafft']), b2b_tsv, aln_fa, 'mafft' ]
        }

    //
    // Join b2b TSV with MMseqs2 pairwise alignment TSV
    //
    ch_b2b_mmseqs = ch_og_b2b
        .join(ch_mmseqs_out)
        .map { meta, b2b_tsv, pairali_tsv ->
            [ meta.plus([mode: 'mmseqs2']), b2b_tsv, pairali_tsv, 'mmseqs2' ]
        }

    //
    // Build distance matrix — runs for both alignment modes on the mixed channel
    //
    B2B_DISTANCE_MATRIX(ch_b2b_mafft.mix(ch_b2b_mmseqs))
    ch_versions = ch_versions.mix(B2B_DISTANCE_MATRIX.out.versions)

    emit:
    distance_matrices = B2B_DISTANCE_MATRIX.out.matrix  // channel: [ val(meta), path(*_b2b_dist.csv) ]
    mafft_alignments  = ch_mafft_out                       // channel: [ val(meta), path(*_aln.fa) ]
    mmseqs_pairali    = ch_mmseqs_out                     // channel: [ val(meta), path(*_pairali.tsv) ]
    versions          = ch_versions
}
