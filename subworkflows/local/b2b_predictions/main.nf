/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B2B_PREDICTIONS subworkflow
    ─────────────────────────────────────────────────────────────────────────────────
    1. Tag every sequence header with its OG ID and pool across all OGs.
    2. Split the pooled FASTA into batches of params.b2b_batch_size sequences;
       each batch becomes one independent Slurm job.
    3. Run b2bTools on every batch in parallel.
    4. Merge all batch TSVs and split back into per-OG files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { POOL_SEQUENCES } from '../../../modules/local/pool_sequences/main'
include { B2BTOOLS       } from '../../../modules/local/b2btools/main'
include { REATTACH_B2B   } from '../../../modules/local/reattach_b2b/main'

workflow B2B_PREDICTIONS {

    take:
    ch_og_fastas    // channel: [ val(meta), path(fasta) ]
    og_dir          // path: directory containing all OG FASTA files

    main:
    ch_versions = Channel.empty()

    //
    // Pool all filtered OG FASTAs (headers renamed to OG_ID|seq_id).
    // Files are read directly from og_dir via glob pattern.
    //
    POOL_SEQUENCES(og_dir)
    ch_versions = ch_versions.mix(POOL_SEQUENCES.out.versions)
    ch_names = POOL_SEQUENCES.out.names_map

    //
    // Split the single pooled FASTA into batches of ~params.b2b_batch_size sequences.
    // splitFasta emits one path per batch file; map adds a synthetic batch meta.
    //
    ch_batches = POOL_SEQUENCES.out.pooled_fasta
        .splitFasta(by: params.b2b_batch_size, file: true)
        .map { batch_file ->
            def batch_id = "batch_${batch_file.baseName}"
            [ [id: batch_id], batch_file ]
        }

    //
    // b2bTools predictions — one Slurm job per batch (fully parallel)
    //
    B2BTOOLS(ch_batches)
    ch_versions = ch_versions.mix(B2BTOOLS.out.versions)

    //
    // Merge all per-batch TSVs into one file, then reattach rows to OGs.
    // keepHeader: true  — copy header from the first TSV
    // skip: 1           — skip header line in all subsequent TSVs
    //
    ch_merged = B2BTOOLS.out.predictions
        .map { meta, tsv -> tsv }
        .collectFile(name: 'all_b2b_predictions.tsv', keepHeader: true, skip: 1)

    REATTACH_B2B(ch_merged, ch_names)
    ch_versions = ch_versions.mix(REATTACH_B2B.out.versions)

    //
    // Reconstruct [ meta, og_b2b_tsv ] tuples from the per-OG output files.
    // File names have the form  <og_id>_b2b.tsv  (written by reattach_b2b.py).
    //
    ch_og_b2b = REATTACH_B2B.out.og_predictions
        .flatten()
        .map { tsv_file ->
            def og_id = tsv_file.baseName.replace('_b2b', '')
            [ [id: og_id], tsv_file ]
        }

    emit:
    og_b2b   = ch_og_b2b    // channel: [ val(meta), path(*_b2b.tsv) ]
    versions = ch_versions
}
