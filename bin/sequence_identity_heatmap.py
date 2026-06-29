#!/usr/bin/env python3
"""Generate heatmap visualization from a pre-computed sequence identity matrix CSV.

This script is separated from the matrix computation to allow independent
matrix generation and visualization processes.

Input:
  - Identity matrix CSV file (from sequence_identity_matrix.py)
  
Output:
  - Heatmap PNG visualization
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate heatmap from sequence identity matrix")
    p.add_argument("--matrix", required=True,
                  help="Input identity matrix CSV file")
    p.add_argument("--og-id", required=True,
                  help="Orthogroup ID for naming output files")
    p.add_argument("--mode", choices=["mafft", "mmseqs2"], required=True,
                  help="Alignment mode: 'mafft' or 'mmseqs2'")
    p.add_argument("--out-heatmap", required=True,
                  help="Output heatmap PNG file path")
    p.add_argument("--cluster-labels", default=None,
                  help="Optional CSV file with cluster labels (sequence_id,cluster columns)")
    p.add_argument("--ref-matrix", default=None,
                  help="Optional reference distance matrix CSV to use for sequence ordering")
    p.add_argument("--sequence-order", default=None,
                  help="Optional file with sequence order (one sequence ID per line)")
    return p.parse_args()


def load_cluster_labels(path: str, seq_ids: list[str]) -> np.ndarray | None:
    """Load cluster labels from CSV file."""
    if path is None:
        return None
    
    try:
        df = pd.read_csv(path, dtype=str)
        if "sequence_id" not in df.columns or "cluster" not in df.columns:
            return None
        
        label_map = dict(zip(df["sequence_id"], df["cluster"]))
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
            file_seq_ids = [line.strip() for line in f if line.strip()]
        
        ordered_seqs = [sid for sid in file_seq_ids if sid in seq_ids]
        
        file_set = set(file_seq_ids)
        for sid in seq_ids:
            if sid not in file_set:
                ordered_seqs.append(sid)
        
        return ordered_seqs
    except Exception:
        return None


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
        sorted_ids = [sid for sid in sequence_order if sid in matrix.index]
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
        sorted_vals = matrix.values
        sorted_ids = matrix.index.tolist()
        sorted_labels = None
    
    n = len(sorted_ids)
    size = max(5, n * 0.35)
    fig, ax = plt.subplots(figsize=(size, size * 0.85))
    
    masked_vals = np.ma.masked_invalid(sorted_vals)
    
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


def main() -> int:
    args = parse_args()
    
    if not Path(args.matrix).exists():
        print(f"Error: Matrix file not found: {args.matrix}", file=sys.stderr)
        return 1
    
    # Load matrix
    try:
        matrix = pd.read_csv(args.matrix, index_col=0)
    except Exception as e:
        print(f"Error reading matrix file: {e}", file=sys.stderr)
        return 1
    
    if matrix.empty:
        print("Error: Matrix is empty", file=sys.stderr)
        return 1
    
    seq_ids = matrix.index.tolist()
    
    # Determine sequence ordering
    sequence_order = None
    
    if args.sequence_order and Path(args.sequence_order).exists():
        sequence_order = load_sequence_order(args.sequence_order, seq_ids)
    elif args.ref_matrix and Path(args.ref_matrix).exists():
        try:
            ref_df = pd.read_csv(args.ref_matrix, index_col=0)
            seq_ids = ref_df.index.tolist()
            if list(ref_df.columns) != seq_ids:
                all_seqs = list(dict.fromkeys(list(ref_df.index) + list(ref_df.columns)))
                seq_ids = sorted(all_seqs)
        except Exception:
            seq_ids = seq_ids
    
    # Load cluster labels
    labels = None
    if not args.sequence_order and not args.ref_matrix:
        labels = load_cluster_labels(args.cluster_labels, seq_ids)
    
    # Generate heatmap
    if sequence_order:
        plot_identity_heatmap(matrix, None, sequence_order, args.og_id, args.mode, args.out_heatmap)
    else:
        plot_identity_heatmap(matrix, labels, None, args.og_id, args.mode, args.out_heatmap)
    
    print(f"Successfully created identity heatmap for {args.og_id} ({args.mode})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
