/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MMSEQS2_EASYSEARCH — all-vs-all pairwise alignment within an OG
    ─────────────────────────────────────────────────────────────────────────────────
    Output columns in pairali.tsv (tab-separated, no header):
      query  target  qaln  taln  qstart  qend  tstart  tend  qlen  tlen  pident  evalue
    Self-hits and direction deduplication are handled in b2b_distance_matrix.py.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MMSEQS2_EASYSEARCH {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.pairali.tsv"), emit: pairali
    path 'versions.yml'                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def og_id = fasta.baseName
    """
    mkdir -p tmp

    mmseqs easy-search \
        "${fasta}" \
        "${fasta}" \
        "${og_id}.pairali.tsv" \
        tmp \
        --format-output "query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,pident,evalue" \
        --threads ${task.cpus} \
        ${args}

    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pairali.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: 15.6f452
    END_VERSIONS
    """
}
