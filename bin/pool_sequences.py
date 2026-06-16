#!/usr/bin/env python3
"""Pool FASTA files and rewrite headers to OG_ID|sequence_id."""
###CONCERN: check if OG>1k seqs, will 1 complete file be created 
from __future__ import annotations

import os
import sys


def parse_fasta(path: str):
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def og_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    if "." in base:
        return base.split(".", 1)[0]
    return base


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: pool_sequences.py <in1.fa> [in2.fa ...] <out.fa>\n")
        return 2

    *input_fastas, out_fa = sys.argv[1:]

    # derive mapping TSV path from output FASTA
    base, _ = os.path.splitext(out_fa)
    map_path = base + ".names.tsv"

    with open(out_fa, "w", encoding="utf-8") as out_handle, open(map_path, "w", encoding="utf-8") as map_handle:
        # header: OG, original sequence name, stripped name (dots removed)
        map_handle.write("OG\toriginal_name\tstripped_name\n")
        for fasta in input_fastas:
            og_id = og_id_from_filename(fasta)
            for seq_id, seq in parse_fasta(fasta):
                clean_id = seq_id.split()[0]
                stripped = clean_id.replace('.', '')
                # write pooled FASTA record
                out_handle.write(f">{og_id}|{clean_id}\n")
                out_handle.write(seq + "\n")
                # write mapping row for later correction step
                map_handle.write(f"{og_id}\t{clean_id}\t{stripped}\n")

    return 0



if __name__ == "__main__":
    raise SystemExit(main())
