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
    chunk_size 

    main:
    ch_versions = Channel.empty()

    // === CRITICAL DEBUG: Check input channels ===
    def count_fastas = 0
    def count_b2b = 0

    ch_og_fastas
        .map { meta, fasta ->
            count_fastas++
            log.info "[ALIGNMENT_DISTANCES] ch_og_fastas: ${meta.id} -> ${fasta} (exists: ${fasta.exists()})"
            [meta, fasta]
        }

    ch_og_b2b
        .map { meta, b2b ->
            count_b2b++
            log.info "[ALIGNMENT_DISTANCES] ch_og_b2b: ${meta.id} -> ${b2b} (exists: ${b2b.exists()})"
            [meta, b2b]
        }

    log.info "[ALIGNMENT_DISTANCES] Received ${count_fastas} FASTA files and ${count_b2b} B2B files"
    //
    // Global MSA — MAFFT with --reorder (run in chunks of params.chunk_size)
    //
    ch_mafft_chunks = ch_og_fastas
        .buffer(chunk_size)
        .map { chunk ->
            def meta = [id: "mafft_chunk_${chunk.hashCode().abs() % 1000000}"]
            def fastas = chunk.collect { it[1] }
            [meta, fastas]
        }// DEBUG: Verify chunks were created
        .map { meta, fastas ->
            log.info "[MAFFT_CHUNKS] Created chunk ${meta.id} with ${fastas.size()} FASTAs"
            [meta, fastas]
        }
    
    MAFFT(ch_mafft_chunks)
    ch_versions = ch_versions.mix(MAFFT.out.versions)

    //
    // Flatten MAFFT outputs back to per-OG for downstream joins
    //
    ch_mafft_out = MAFFT.out.alignment
        .flatMap { meta, aln_files ->
            // aln_files is a list of path objects
            def files = aln_files instanceof List ? aln_files : [aln_files]
            files.collect { aln_file ->
                def og_id = aln_file.baseName.replace('.aln.fa', '')
                [ [id: og_id], aln_file ]
            }
        }

    //
    // Pairwise alignment — MMseqs2 easy-search all-vs-all per OG (run in chunks)
    //
    ch_mmseqs_chunks = ch_og_fastas
        .buffer(chunk_size)
        .map { chunk ->
            def meta = [id: "mmseqs_chunk_${chunk.hashCode().abs() % 1000000}"]
            def fastas = chunk.collect { it[1] }
            [meta, fastas]
        }// DEBUG: Verify chunks were created
        .map { meta, fastas ->
            log.info "[MMSEQS_CHUNKS] Created chunk ${meta.id} with ${fastas.size()} FASTAs"
            [meta, fastas]
        }
    
    MMSEQS2_EASYSEARCH(ch_mmseqs_chunks)
    ch_versions = ch_versions.mix(MMSEQS2_EASYSEARCH.out.versions)

    //
    // Flatten MMSEQS2 outputs back to per-OG for downstream joins
    //
    ch_mmseqs_out = MMSEQS2_EASYSEARCH.out.pairali
        .flatMap { meta, pairali_files ->
            def files = pairali_files instanceof List ? pairali_files : [pairali_files]
            files.collect { pairali_file ->
                def og_id = pairali_file.baseName.replace('.pairali.tsv', '')
                [ [id: og_id], pairali_file ]
            }
        }

    //
    // Join b2b TSV with MAFFT alignment (key: meta.id)
    // Augment meta with mode so output filenames and PCA labels are unambiguous.
    //
    ch_b2b_mafft = ch_og_b2b
        .join(ch_mafft_out)
        .map { meta, b2b_tsv, aln_fa ->
            [ meta + [mode: 'mafft'], b2b_tsv, aln_fa, 'mafft' ]
        }

    //
    // Join b2b TSV with MMseqs2 pairwise alignment TSV
    //
    ch_b2b_mmseqs = ch_og_b2b
        .join(ch_mmseqs_out)
        .map { meta, b2b_tsv, pairali_tsv ->
            [ meta + [mode: 'mmseqs2'], b2b_tsv, pairali_tsv, 'mmseqs2' ]
        }

    //
    // Build distance matrix — runs for both alignment modes on the mixed channel
    //
    B2B_DISTANCE_MATRIX(ch_b2b_mafft.mix(ch_b2b_mmseqs))
    ch_versions = ch_versions.mix(B2B_DISTANCE_MATRIX.out.versions)

    emit:
    distance_matrices = B2B_DISTANCE_MATRIX.out.matrix  // channel: [ val(meta), path(*_b2b_dist.csv) ]
    versions          = ch_versions
}
