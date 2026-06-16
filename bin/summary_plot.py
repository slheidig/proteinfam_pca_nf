#!/usr/bin/env python3
"""Aggregate per-OG clustering outputs and generate a summary bar plot."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--clusters", nargs="+", required=True, help="One or more *_clusters.csv files")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-hist-csv", required=True)
    p.add_argument("--out-plot", required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    frames = []
    for fp in args.clusters:
        path = Path(fp)
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path)
        if {"og_id", "mode", "cluster"}.issubset(df.columns):
            frames.append(df)

    if not frames:
        out_df = pd.DataFrame(columns=["og_id", "mode", "n_clusters"])
        out_df.to_csv(args.out_csv, index=False)
        pd.DataFrame(columns=["mode", "cluster_bin", "n_ogs"]).to_csv(args.out_hist_csv, index=False)

        plt.figure(figsize=(8, 4))
        plt.title("No cluster files found")
        plt.tight_layout()
        plt.savefig(args.out_plot, dpi=150)
        return 0

    all_df = pd.concat(frames, ignore_index=True)
    summary = (
        all_df.groupby(["og_id", "mode"], as_index=False)["cluster"]
        .nunique()
        .rename(columns={"cluster": "n_clusters"})
        .sort_values(["mode", "n_clusters", "og_id"], ascending=[True, False, True])
    )
    # Round numeric columns to 4 decimals for CSV outputs
    summary = summary.copy()
    summary = summary.round(4)
    summary.to_csv(args.out_csv, index=False)

    modes = summary["mode"].drop_duplicates().tolist()
    # Exclude the '1' cluster bin from summary histograms/CSV (show only >=2) because thats where k-means starts scan
    bin_labels = [str(i) for i in range(2, 11)] + ["10+"]

    hist_rows = []
    for mode in modes:
        sub = summary[summary["mode"] == mode].copy()
        sub["cluster_bin"] = sub["n_clusters"].apply(lambda x: "10+" if int(x) > 10 else str(int(x)))
        counts = sub["cluster_bin"].value_counts()
        for b in bin_labels:
            hist_rows.append({"mode": mode, "cluster_bin": b, "n_ogs": int(counts.get(b, 0))})
    hist_df = pd.DataFrame(hist_rows)
    hist_df = hist_df.copy()
    hist_df = hist_df.round(4)
    hist_df.to_csv(args.out_hist_csv, index=False)

    n_modes = max(len(modes), 1)

    fig, axes = plt.subplots(n_modes, 1, figsize=(10, 4 * n_modes), sharey=True)
    if n_modes == 1:
        axes = [axes]

    for ax, mode in zip(axes, modes):
        sub_hist = hist_df[hist_df["mode"] == mode].copy()
        sub_hist["cluster_bin"] = pd.Categorical(sub_hist["cluster_bin"], categories=bin_labels, ordered=True)
        sub_hist = sub_hist.sort_values("cluster_bin")
        ax.bar(sub_hist["cluster_bin"].astype(str), sub_hist["n_ogs"], width=0.85)
        ax.set_title(f"Histogram of internal clusters ({mode})")
        ax.set_xlabel("Internal clusters per OG (>=2)")
        ax.set_ylabel("Number of OGs")

    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=150)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
