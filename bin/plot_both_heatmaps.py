#!/usr/bin/env python3
"""Plot both b2b distance and sequence identity heatmaps using the same sequence order.

This script is independent from PCA calculation. It takes pre-computed matrices,
cluster labels, and sequence order, and generates both heatmaps with consistent
sorting.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--b2b-matrix", required=True, help="B2B distance matrix CSV file")
    p.add_argument("--seq-id-matrix", required=True, help="Sequence identity matrix CSV file")
    p.add_argument("--cluster-labels", required=True, help="Cluster labels CSV file")
    p.add_argument("--sequence-order", required=True, help="Sequence order file (one ID per line)")
    p.add_argument("--og-id", required=True, help="Orthogroup ID")
    p.add_argument("--mode", required=True, help="Alignment mode (mafft/mmseqs2)")
    p.add_argument("--out-b2b-heatmap", required=True, help="Output path for b2b distance heatmap")
    p.add_argument("--out-seq-heatmap", required=True, help="Output path for sequence identity heatmap")
    return p.parse_args()


def load_sequence_order(path: str, all_seq_ids: list[str]) -> list[str]:
    """Load sequence order from file, filtering to available sequences."""
    try:
        with open(path, "r") as f:
            ordered = [line.strip() for line in f if line.strip()]
        # Filter to sequences present in all_seq_ids, preserving order
        return [s for s in ordered if s in all_seq_ids]
    except Exception:
        return all_seq_ids


def load_cluster_labels(path: str, seq_ids: list[str]) -> np.ndarray:
    """Load cluster labels from CSV file."""
    try:
        df = pd.read_csv(path, dtype=str)
        if "sequence_id" not in df.columns or "cluster" not in df.columns:
            return np.zeros(len(seq_ids), dtype=int)
        label_map = dict(zip(df["sequence_id"], df["cluster"]))
        return np.array([label_map.get(sid, 0) for sid in seq_ids], dtype=int)
    except Exception:
        return np.zeros(len(seq_ids), dtype=int)


def load_matrix(path: str) -> pd.DataFrame:
    """Load a matrix CSV file."""
    return pd.read_csv(path, index_col=0)


def plot_single_heatmap(
    mat: pd.DataFrame,
    seq_order: list[str],
    og_id: str,
    mode: str,
    out_path: str,
    matrix_type: str,
    labels: np.ndarray | None = None,
) -> None:
    """Plot a single heatmap with the given sequence order.
    
    Args:
        matrix_type: 'b2b_distance' (0-0.15 scale) or 'seq_identity' (0-100 scale)
        labels: Optional cluster labels for drawing cluster boundaries
    """
    # Reorder matrix by sequence order
    # Filter seq_order to only include sequences in matrix
    valid_order = [s for s in seq_order if s in mat.index]
    # Add any missing sequences at the end
    for sid in mat.index:
        if sid not in valid_order:
            valid_order.append(sid)
    
    sorted_mat = mat.reindex(index=valid_order, columns=valid_order)
    sorted_vals = sorted_mat.values
    
    # If labels are provided, filter and reorder them to match valid_order
    sorted_labels = None
    if labels is not None and len(labels) == len(mat.index):
        label_map = dict(zip(mat.index, labels))
        sorted_labels = np.array([label_map.get(sid, 0) for sid in valid_order], dtype=int)
    
    # Mask NaN values
    masked_vals = np.ma.masked_invalid(sorted_vals)
    
    n = len(valid_order)
    size = max(5, n * 0.35)
    fig, ax = plt.subplots(figsize=(size, size * 0.85))

    cmap = plt.colormaps["viridis_r"].copy()
    cmap.set_bad(color="lightgrey")
    
    # Set vmin/vmax and colorbar label based on matrix type
    if matrix_type == "seq_identity":
        im = ax.imshow(masked_vals, aspect="auto", cmap=cmap, interpolation="none", vmin=0, vmax=100)
        cbar_label = "Sequence identity (%)"
    else:  # b2b_distance
        im = ax.imshow(masked_vals, aspect="auto", cmap=cmap, interpolation="none", vmin=0, vmax=0.15)
        cbar_label = "B2B distance"
    
    plt.colorbar(im, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)

    ax.set_xticks(range(n))
    ax.set_xticklabels(valid_order, rotation=90, fontsize=4)
    ax.set_yticks(range(n))
    ax.set_yticklabels(valid_order, fontsize=4)
    ax.tick_params(axis="x", pad=0.1, length=0)
    ax.tick_params(axis="y", pad=0.1, length=0)

    # Draw cluster boundary lines if labels are provided
    if sorted_labels is not None and len(sorted_labels) > 0:
        prev = sorted_labels[0]
        for i, lab in enumerate(sorted_labels[1:], start=1):
            if lab != prev:
                ax.axhline(i - 0.5, color="red", lw=1.0)
                ax.axvline(i - 0.5, color="red", lw=1.0)
                prev = lab

    # Determine matrix name for title
    matrix_name = "seq identity" if matrix_type == "seq_identity" else "b2b distance"
    ax.set_title(f"{og_id} | {mode} | {matrix_name} matrix (sorted by PCA clusters)")
    
    plt.tight_layout()
    size = max(5, min(n * 0.35, 20))
    plt.savefig(out_path, dpi=250)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    
    # Ensure output directory exists
    for path in [args.out_b2b_heatmap, args.out_seq_heatmap]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    
    # Load inputs
    b2b_mat = load_matrix(args.b2b_matrix)
    seq_mat = load_matrix(args.seq_id_matrix)
    seq_order = load_sequence_order(args.sequence_order, b2b_mat.index.tolist())
    cluster_labels = load_cluster_labels(args.cluster_labels, b2b_mat.index.tolist())
    
    # Plot both heatmaps using the same sequence order and cluster labels
    plot_single_heatmap(
        b2b_mat,
        seq_order,
        args.og_id,
        args.mode,
        args.out_b2b_heatmap,
        matrix_type="b2b_distance",
        labels=cluster_labels,
    )
    
    plot_single_heatmap(
        seq_mat,
        seq_order,
        args.og_id,
        args.mode,
        args.out_seq_heatmap,
        matrix_type="seq_identity",
        labels=cluster_labels,
    )
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
