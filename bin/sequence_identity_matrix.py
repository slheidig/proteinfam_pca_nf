#!/usr/bin/env python3
"""Compute pairwise sequence identity matrices from MAFFT MSA or MMseqs2 pairwise alignments.

For MAFFT: Calculate sequence identity from the MSA for every pair of sequences.
  - Skip full gap positions (-,-)
  - Gap positions count as mismatch (0 points)
  - Each match = 1 point
  - Sum points and divide by number of compared positions
  - Scale to 0-100%

For MMseqs2: Read pairwise sequence identity from the pident column (second-last column).
  - If a pair is missing, leave as NaN
  - NaN values are masked grey in the heatmap

Output:
  - Full NxN matrix CSV: {og_id}_{mode}_identity_matrix.csv
  - Heatmap PNG: {og_id}_{mode}_identity_heatmap.png
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute pairwise sequence identity matrices")
    p.add_argument("--ali", required=True, 
                  help="Alignment input file (MAFFT fasta or MMseqs2 pairali TSV)")
    p.add_argument("--mode", choices=["mafft", "mmseqs2"], required=True,
                  help="Alignment mode: 'mafft' for MSA FASTA, 'mmseqs2' for pairali TSV")
    p.add_argument("--og-id", required=True,
                  help="Orthogroup ID for naming output files")
    p.add_argument("--out-matrix", required=True,
                  help="Output CSV matrix file path")
    p.add_argument("--out-heatmap", required=True,
                  help="Output heatmap PNG file path")
    p.add_argument("--cluster-labels", default=None,
                  help="Optional CSV file with cluster labels (sequence_id,cluster columns)")
    p.add_argument("--ref-matrix", default=None,
                  help="Optional reference distance matrix CSV to use for sequence ordering")
    p.add_argument("--sequence-order", default=None,
                  help="Optional file with sequence order (one sequence ID per line) to use for heatmap sorting")
    return p.parse_args()


def parse_fasta_alignment(path: str) -> dict[str, str]:
    """Parse a FASTA alignment file and return dict of sequence_id -> aligned_sequence."""
    ali = {}
    sid = None
    parts = []
    
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


def parse_mmseqs_pairali(path: str) -> dict[frozenset, float]:
    """Parse MMseqs2 pairali TSV and return dict of frozenset({a,b}) -> pident percentage.
    
    Expected columns (tab-separated, no header):
    query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,pident,evalue
    
    Returns only non-self pairs, with pident as float percentage (0-100).
    """
    try:
        df = pd.read_csv(path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return {}
    
    # Column indices: pident is at index 10 (11th column)
    if df.shape[1] <= 10:
        return {}
    
    # Filter out self-hits and create canonical unordered pairs
    mask = df.iloc[:, 0] != df.iloc[:, 1]  # query != target
    pair_df = df.loc[mask, [0, 1, 10]].copy()
    pair_df.columns = ["query", "target", "pident"]
    
    pair_map = {}
    for row in pair_df.itertuples(index=False):
        key = frozenset((row.query, row.target))
        try:
            pident = float(row.pident)
        except (ValueError, TypeError):
            continue
        # Keep first occurrence (or could average; using first for simplicity)
        if key not in pair_map:
            pair_map[key] = pident
    
    return pair_map


def calculate_sequence_identity_from_msa(
    ali: dict[str, str],
    seq_ids: list[str]
) -> dict[tuple[str, str], float]:
    """Calculate pairwise sequence identity from MSA alignment.
    
    For each pair (a, b):
      - Skip positions where both are gaps (-,-)
      - Gap in one sequence counts as mismatch (0 points)
      - Match = 1 point
      - identity = (sum of matches / number of compared positions) * 100
    
    Returns dict keyed by (a, b) with a < b for uniqueness.
    """
    pairs = list(itertools.combinations(seq_ids, 2))
    results = {}
    
    for a, b in pairs:
        if a not in ali or b not in ali:
            results[(a, b)] = float('nan')
            continue
        
        seq_a = ali[a]
        seq_b = ali[b]
        
        if len(seq_a) != len(seq_b):
            results[(a, b)] = float('nan')
            continue
        
        matches = 0
        compared = 0
        
        for ca, cb in zip(seq_a, seq_b):
            a_gap = ca == "-"
            b_gap = cb == "-"
            
            # Skip full gap positions
            if a_gap and b_gap:
                continue
            
            # Count position as compared
            compared += 1
            
            # Match = both are same non-gap character
            if not a_gap and not b_gap and ca == cb:
                matches += 1
        
        if compared == 0:
            results[(a, b)] = float('nan')
        else:
            results[(a, b)] = (matches / compared) * 100.0
    
    return results


def build_matrix(
    pair_results: dict[tuple[str, str], float] | dict[frozenset, float],
    seq_ids: list[str],
    mode: str
) -> pd.DataFrame:
    """Build a full NxN matrix DataFrame from pairwise results."""
    n = len(seq_ids)
    matrix = pd.DataFrame(np.nan, index=seq_ids, columns=seq_ids, dtype=float)
    
    # Set diagonal to 100% (identity with self)
    np.fill_diagonal(matrix.values, 100.0)
    
    for key, value in pair_results.items():
        if mode == "mafft":
            a, b = key
        else:  # mmseqs2
            items = list(key)
            a, b = items[0], items[1]
        
        matrix.loc[a, b] = value
        matrix.loc[b, a] = value
    
    return matrix


def plot_identity_heatmap(
    matrix: pd.DataFrame,
    labels: np.ndarray | None,
    sequence_order: list[str] | None,
    og_id: str,
    mode: str,
    out_path: str,
) -> None:
    """Save a heatmap of the identity matrix with rows/cols sorted by cluster label or custom order.
    
    NaN values are masked and shown in grey.
    
    Args:
        matrix: The identity matrix DataFrame
        labels: Optional cluster labels for sorting
        sequence_order: Optional explicit sequence order (overrides labels)
        og_id: Orthogroup ID
        mode: Alignment mode
        out_path: Output file path for the heatmap
    """
    if sequence_order is not None:
        # Use the explicit sequence order
        # Filter to only include sequences that are in the matrix
        sorted_ids = [sid for sid in sequence_order if sid in matrix.index]
        # Add any missing sequences from matrix at the end
        for sid in matrix.index:
            if sid not in sorted_ids:
                sorted_ids.append(sid)
        
        sorted_vals = matrix.reindex(index=sorted_ids, columns=sorted_ids).values
        sorted_labels = None
    elif labels is not None:
        order = np.argsort(labels, kind="stable")
        sorted_vals = matrix.values[np.ix_(order, order)]
        sorted_ids = [matrix.index[i] for i in order]
        sorted_labels = labels[order]
    else:
        # No cluster labels, use original order
        sorted_vals = matrix.values
        sorted_ids = matrix.index.tolist()
        sorted_labels = None
    
    n = len(sorted_ids)
    size = max(5, n * 0.35)
    fig, ax = plt.subplots(figsize=(size, size * 0.85))
    
    # Create masked array for NaN values
    masked_vals = np.ma.masked_invalid(sorted_vals)
    
    # Use viridis colormap (reversed for identity: high identity = dark)
    cmap = plt.cm.get_cmap("viridis_r")
    try:
        cmap = cmap.copy()
    except Exception:
        pass
    cmap.set_bad(color="lightgrey")
    
    im = ax.imshow(masked_vals, aspect="auto", cmap=cmap, 
                  interpolation="none", vmin=0, vmax=100)
    plt.colorbar(im, ax=ax, label="Sequence identity (%)", fraction=0.046, pad=0.04)
    
    ax.set_xticks(range(n))
    ax.set_xticklabels(sorted_ids, rotation=90, fontsize=4)
    ax.set_yticks(range(n))
    ax.set_yticklabels(sorted_ids, fontsize=4)
    ax.tick_params(axis="x", pad=0.1, length=0)
    ax.tick_params(axis="y", pad=0.1, length=0)
    
    # Red lines at cluster boundaries if labels are provided
    if sorted_labels is not None:
        prev = sorted_labels[0]
        for i, lab in enumerate(sorted_labels[1:], start=1):
            if lab != prev:
                ax.axhline(i - 0.5, color="red", lw=1.0)
                ax.axvline(i - 0.5, color="red", lw=1.0)
                prev = lab
    
    ax.set_title(f"{og_id} | {mode} | identity matrix")
    plt.tight_layout()
    plt.savefig(out_path, dpi=350)
    plt.close(fig)


def load_cluster_labels(path: str, seq_ids: list[str]) -> np.ndarray | None:
    """Load cluster labels from CSV file."""
    if path is None:
        return None
    
    try:
        df = pd.read_csv(path, dtype=str)
        if "sequence_id" not in df.columns or "cluster" not in df.columns:
            return None
        
        # Create mapping from sequence_id to cluster
        label_map = dict(zip(df["sequence_id"], df["cluster"]))
        
        # Map seq_ids to labels, filling missing with -1
        labels = []
        for sid in seq_ids:
            labels.append(label_map.get(sid, -1))
        
        return np.array(labels, dtype=int)
    except Exception:
        return None


def load_sequence_order(path: str, seq_ids: list[str]) -> list[str] | None:
    """Load sequence order from file (one sequence ID per line).
    
    Returns a list of sequence IDs in the desired order.
    Only includes sequences that are present in seq_ids.
    """
    if path is None:
        return None
    
    try:
        with open(path, "r") as f:
            # Read all sequence IDs from the file, stripping whitespace
            file_seq_ids = [line.strip() for line in f if line.strip()]
        
        # Filter to only include sequences that are in our seq_ids list
        # and preserve the order from the file
        ordered_seqs = [sid for sid in file_seq_ids if sid in seq_ids]
        
        # Add any missing sequences at the end (maintaining stability)
        file_set = set(file_seq_ids)
        for sid in seq_ids:
            if sid not in file_set:
                ordered_seqs.append(sid)
        
        return ordered_seqs
    except Exception:
        return None


def main() -> int:
    args = parse_args()
    
    # Validate inputs
    if not Path(args.ali).exists():
        print(f"Error: Alignment file not found: {args.ali}", file=sys.stderr)
        return 1
    
    # Read alignment data
    if args.mode == "mafft":
        ali = parse_fasta_alignment(args.ali)
        seq_ids_from_ali = sorted(ali.keys())
        
        if not seq_ids_from_ali:
            print("Error: No sequences found in MAFFT alignment", file=sys.stderr)
            return 1
        
        # Calculate pairwise identity from MSA
        pair_results = calculate_sequence_identity_from_msa(ali, seq_ids_from_ali)
        
    else:  # mmseqs2
        pair_map = parse_mmseqs_pairali(args.ali)
        
        # Extract all unique sequence IDs from the pairali file
        all_ids = set()
        try:
            df = pd.read_csv(args.ali, sep="\t", header=None)
            if df.shape[1] >= 2:
                all_ids = set(df.iloc[:, 0].dropna().unique()) | set(df.iloc[:, 1].dropna().unique())
        except Exception:
            all_ids = set()
        
        if not all_ids:
            # Try to get seq_ids from cluster labels if available
            if args.cluster_labels:
                try:
                    df = pd.read_csv(args.cluster_labels, dtype=str)
                    if "sequence_id" in df.columns:
                        all_ids = set(df["sequence_id"].dropna().unique())
                except Exception:
                    pass
        
        if not all_ids:
            print("Error: No sequence IDs found in MMseqs2 output", file=sys.stderr)
            return 1
        
        seq_ids_from_ali = sorted(all_ids)
        pair_results = pair_map  # Already in frozenset format
    
    # Determine sequence ordering
    # Priority: 1) explicit sequence order file, 2) reference matrix, 3) cluster labels, 4) sorted from alignment
    sequence_order = None
    
    if args.sequence_order and Path(args.sequence_order).exists():
        # Load explicit sequence order from file
        sequence_order = load_sequence_order(args.sequence_order, seq_ids_from_ali)
        if sequence_order:
            seq_ids = sequence_order
        else:
            seq_ids = seq_ids_from_ali
    elif args.ref_matrix and Path(args.ref_matrix).exists():
        # Use the row/column order from the reference distance matrix
        try:
            ref_df = pd.read_csv(args.ref_matrix, index_col=0)
            seq_ids = ref_df.index.tolist()
            # Also check columns match index
            if list(ref_df.columns) != seq_ids:
                # If columns differ, use sorted union
                all_seqs = list(dict.fromkeys(list(ref_df.index) + list(ref_df.columns)))
                seq_ids = sorted(all_seqs)
        except Exception:
            # Fall back to sorted from alignment
            seq_ids = seq_ids_from_ali
    else:
        seq_ids = seq_ids_from_ali
    
    # Filter pair_results to only include sequences in our final seq_ids list
    if args.mode == "mafft":
        # pair_results is dict of (a, b) -> value
        filtered_results = {k: v for k, v in pair_results.items() 
                          if k[0] in seq_ids and k[1] in seq_ids}
    else:
        # pair_results is dict of frozenset({a, b}) -> value
        filtered_results = {k: v for k, v in pair_results.items() 
                          if all(s in seq_ids for s in k)}
    
    # Build full matrix with the determined sequence order
    matrix = build_matrix(filtered_results, seq_ids, args.mode)
    
    # Load cluster labels for sorting (if no sequence_order or ref_matrix was used)
    labels = None
    if not args.sequence_order and not args.ref_matrix:
        labels = load_cluster_labels(args.cluster_labels, seq_ids)
    
    # Save matrix CSV

    matrix = matrix.round(4)
    matrix.to_csv(args.out_matrix, index=True)
    
    # Generate heatmap
    # If we have an explicit sequence_order, use it; otherwise use cluster labels
    if sequence_order:
        plot_identity_heatmap(matrix, None, sequence_order, args.og_id, args.mode, args.out_heatmap)
    else:
        plot_identity_heatmap(matrix, labels, None, args.og_id, args.mode, args.out_heatmap)
    
    print(f"Successfully created identity matrix and heatmap for {args.og_id} ({args.mode})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
