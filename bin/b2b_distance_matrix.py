#!/usr/bin/env python3
"""Compute pairwise b2b backbone distance matrices from MAFFT or MMseqs2 alignments."""

from __future__ import annotations

import argparse
import itertools
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--b2b", required=True, help="Per-OG b2b TSV")
    p.add_argument("--ali", required=True, help="Alignment input file (MAFFT fasta or MMseqs2 pairali TSV)")
    p.add_argument("--mode", choices=["mafft", "mmseqs2"], required=True)
    p.add_argument("--og-id", required=True)
    p.add_argument("--out-matrix", required=True)
    p.add_argument("--out-meta", required=True)
    return p.parse_args()


def read_b2b_backbone(path: str) -> Dict[str, List[float]]:
    df = pd.read_csv(path, sep="\t")
    required = {"sequence_id", "residue_index", "backbone"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required b2b columns: {sorted(missing)}")

    # Group once in pandas, then build dense backbone vectors per sequence.
    backbone_by_seq: Dict[str, List[float]] = {}
    for sid, group in df.groupby("sequence_id", sort=False):
        if group.empty:
            backbone_by_seq[str(sid)] = []
            continue
        # Filter out NaN and infinite values from residue_index
        group = group.dropna(subset=["residue_index"])
        group = group[~group["residue_index"].isin([float('inf'), float('-inf')])]
        if group.empty:
            backbone_by_seq[str(sid)] = []
            continue
        positions = group["residue_index"].astype(int).to_numpy()
        values = group["backbone"].astype(float).to_numpy()
        max_pos = int(positions.max())
        arr = np.full(max_pos, np.nan, dtype=float)
        arr[positions - 1] = values
        backbone_by_seq[str(sid)] = arr.tolist()
    return backbone_by_seq


def parse_fasta_alignment(path: str) -> Dict[str, str]:
    ali: Dict[str, str] = {}
    sid = None
    parts: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    ali[sid] = "".join(parts)
                sid = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if sid is not None:
            ali[sid] = "".join(parts)
    return ali


def parse_mmseqs_pairali(path: str) -> Dict[frozenset, Tuple[str, str, str, str]]:
    """Return canonical pair alignments keyed by frozenset({a,b}).

    Keep first seen row per unordered pair, drop self-hits.
    Expected columns (no header):
    query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,pident,evalue
    """
    pair_df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["query", "target", "qaln", "taln"],
    )
    pair_df = pair_df[pair_df["query"] != pair_df["target"]].drop_duplicates(
        subset=["query", "target"],
        keep="first",
    )

    pair_map: Dict[frozenset, Tuple[str, str, str, str]] = {}
    for row in pair_df.itertuples(index=False):
        key = frozenset((row.query, row.target))
        if key not in pair_map:
            pair_map[key] = (row.query, row.target, row.qaln, row.taln)
    return pair_map


def distance_from_aligned_strings(
    ali_a: str,
    ali_b: str,
    b2b_a: List[float],
    b2b_b: List[float],
) -> Tuple[float, int, int, int, int]:
    """Map b2b values onto alignment columns and compute normalized mean abs diff.

    Returns (distance, seq_a_len, seq_b_len, ali_len, n_compared).
    """
    if len(ali_a) != len(ali_b):
        raise ValueError("Aligned strings have unequal lengths")

    seq_a_len = sum(1 for c in ali_a if c != "-")
    seq_b_len = sum(1 for c in ali_b if c != "-")
    ali_len = len(ali_a)

    ai = 0
    bi = 0
    total = 0.0
    n = 0

    for ca, cb in zip(ali_a, ali_b):
        a_gap = ca == "-"
        b_gap = cb == "-"

        a_val = None
        b_val = None

        if not a_gap:
            if ai < len(b2b_a):
                a_val = b2b_a[ai]
            ai += 1
        if not b_gap:
            if bi < len(b2b_b):
                b_val = b2b_b[bi]
            bi += 1

        if a_gap or b_gap:
            continue
        if a_val is None or b_val is None:
            continue
        if np.isnan(a_val) or np.isnan(b_val):
            continue

        total += abs(a_val - b_val)
        n += 1

    if n == 0:
        return float("nan"), seq_a_len, seq_b_len, ali_len, 0
    return total / n, seq_a_len, seq_b_len, ali_len, n


def main() -> int:
    args = parse_args()

    backbone = read_b2b_backbone(args.b2b)
    seq_ids = sorted(backbone.keys())

    # Pair list: unordered, no self-pairs, no duplicate direction.
    pairs = list(itertools.combinations(seq_ids, 2))

    matrix = pd.DataFrame(np.nan, index=seq_ids, columns=seq_ids, dtype=float)
    np.fill_diagonal(matrix.values, 0.0)

    meta_rows = []

    if args.mode == "mafft":
        ali = parse_fasta_alignment(args.ali)

        for a, b in pairs:
            if a not in ali or b not in ali:
                continue
            dist, la, lb, lali, ncomp = distance_from_aligned_strings(
                ali[a], ali[b], backbone[a], backbone[b]
            )
            matrix.loc[a, b] = dist
            matrix.loc[b, a] = dist
            meta_rows.append(
                {
                    "og_id": args.og_id,
                    "mode": args.mode,
                    "seq_a": a,
                    "seq_b": b,
                    "seq_a_len": la,
                    "seq_b_len": lb,
                    "ali_len": lali,
                    "compared_residues": ncomp,
                    "distance": dist,
                }
            )

    else:
        pairali = parse_mmseqs_pairali(args.ali)

        for a, b in pairs:
            key = frozenset((a, b))
            if key not in pairali:
                continue
            q, t, qaln, taln = pairali[key]

            # pairali row can come as a->b or b->a; reorder to match (a,b)
            if q == a and t == b:
                ali_a, ali_b = qaln, taln
            elif q == b and t == a:
                ali_a, ali_b = taln, qaln
            else:
                # Defensive fallback for unusual IDs
                continue

            dist, la, lb, lali, ncomp = distance_from_aligned_strings(
                ali_a, ali_b, backbone[a], backbone[b]
            )
            matrix.loc[a, b] = dist
            matrix.loc[b, a] = dist
            meta_rows.append(
                {
                    "og_id": args.og_id,
                    "mode": args.mode,
                    "seq_a": a,
                    "seq_b": b,
                    "seq_a_len": la,
                    "seq_b_len": lb,
                    "ali_len": lali,
                    "compared_residues": ncomp,
                    "distance": dist,
                }
            )

    # Round matrix values to 4 decimals for consistent CSV outputs
    matrix.round(4).to_csv(args.out_matrix, index=True)

    meta_df = pd.DataFrame(meta_rows)
    if meta_df.empty:
        meta_df = pd.DataFrame(
            columns=[
                "og_id",
                "mode",
                "seq_a",
                "seq_b",
                "seq_a_len",
                "seq_b_len",
                "ali_len",
                "compared_residues",
                "distance",
            ]
        )
    # Round numeric columns to 4 decimals before writing meta TSV
    meta_df = meta_df.copy()
    meta_df = meta_df.round(4)
    meta_df.to_csv(args.out_meta, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
