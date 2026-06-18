#!/usr/bin/env python3
"""PCA and silhouette-based clustering for one OG distance matrix."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--og-id", required=True)
    p.add_argument("--mode", required=True)
    p.add_argument("--out-plot", required=True)
    p.add_argument("--out-clusters", required=True)
    p.add_argument("--out-meta", required=True)
    p.add_argument("--out-heatmap", required=True)
    p.add_argument("--out-sequence-order", required=True,
                   help="Output file path for sequence order used in heatmap")
    p.add_argument("--external-labels", default=None,
                   help="TSV file with columns og_id,gene_id,<label1>,<label2>,... for external labels")
    p.add_argument("--summary-dir", default="test_results/summary",
                   help="Directory to append per-OG NaN summary CSV (default: test_results/summary)")
    p.add_argument("--n-perm", type=int, default=999)
    return p.parse_args()


def silhouette_pvalue(x: np.ndarray, labels: np.ndarray, observed_sil: float, n_perm: int = 999) -> float:
    """Permutation test: fraction of random label assignments with silhouette >= observed.

    Null hypothesis: the cluster labels are no better than a random grouping of the
    same size.  Returns (count + 1) / (n_perm + 1) so the p-value is never exactly 0.
    """
    if len(set(labels)) < 2 or np.isnan(observed_sil):
        return float("nan")
    rng = np.random.default_rng(42)
    count = 0
    for _ in range(n_perm):
        perm_labels = rng.permutation(labels)
        if len(set(perm_labels)) < 2:
            continue
        s = silhouette_score(x, perm_labels)
        if s >= observed_sil:
            count += 1
    return (count + 1) / (n_perm + 1)


def plot_heatmap(
    mat: pd.DataFrame,
    labels: np.ndarray,
    og_id: str,
    mode: str,
    out_path: str,
    out_sequence_order_path: str | None = None,
) -> None:
    """Save a heatmap of the distance matrix with rows/cols sorted by cluster label.
    
    If out_sequence_order_path is provided, also saves the sequence order to that file.
    """
    order = np.argsort(labels, kind="stable")
    sorted_vals = mat.values[np.ix_(order, order)]
    # Mask NaN values so they are rendered in a distinct color (grey)
    masked_vals = np.ma.masked_invalid(sorted_vals)
    sorted_ids = [mat.index[i] for i in order]
    sorted_labels = labels[order]
    
    # Save sequence order to file if requested
    if out_sequence_order_path is not None:
        with open(out_sequence_order_path, "w") as f:
            for sid in sorted_ids:
                f.write(f"{sid}\n")

    n = len(sorted_ids)
    size = max(5, n * 0.35)
    fig, ax = plt.subplots(figsize=(size, size * 0.85))

    cmap = plt.cm.get_cmap("viridis_r")
    try:
        cmap = cmap.copy()
    except Exception:
        # Some colormap objects may not support copy; fall back to mutating
        pass
    cmap.set_bad(color="lightgrey")
    im = ax.imshow(masked_vals, aspect="auto", cmap=cmap, interpolation="none")
    plt.colorbar(im, ax=ax, label="B2B distance", fraction=0.046, pad=0.04)

    ax.set_xticks(range(n))
    ax.set_xticklabels(sorted_ids, rotation=90, fontsize=4)
    ax.set_yticks(range(n))
    ax.set_yticklabels(sorted_ids, fontsize=4)
    ax.tick_params(axis="x", pad=0.1, length=0)
    ax.tick_params(axis="y", pad=0.1, length=0)

    # Red lines at cluster boundaries
    prev = sorted_labels[0]
    for i, lab in enumerate(sorted_labels[1:], start=1):
        if lab != prev:
            ax.axhline(i - 0.5, color="red", lw=1.0)
            ax.axvline(i - 0.5, color="red", lw=1.0)
            prev = lab

    ax.set_title(f"{og_id} | {mode} | distance matrix (sorted by cluster)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=350)
    plt.close(fig)


def map_categories(labels: pd.Series) -> np.ndarray:
    """Map arbitrary categorical labels to integers for coloring."""
    categories = list(OrderedDict.fromkeys(labels.astype(str)))
    mapping = {cat: idx for idx, cat in enumerate(categories)}
    return labels.map(mapping).to_numpy(dtype=int)


def map_categories_with_names(labels: pd.Series) -> tuple[np.ndarray, list[str]]:
    """Map categorical labels to integers and also return category names in order.

    Returns (numeric_array, categories_list).
    """
    categories = list(OrderedDict.fromkeys(labels.astype(str)))
    mapping = {cat: idx for idx, cat in enumerate(categories)}
    return labels.map(mapping).to_numpy(dtype=int), categories


def find_sequence_genome(seq_ids: list[str], genomes: pd.Index) -> pd.Series:
    """Resolve each sequence ID to a genome label by longest prefix match."""
    sorted_genomes = sorted((g for g in genomes if pd.notna(g)), key=len, reverse=True)
    seq_series = pd.Series(seq_ids, dtype="string", name="sequence_id")
    matched = pd.Series([pd.NA] * len(seq_series), dtype="string")

    # Longest-prefix wins to avoid partial matches (e.g. Run1_Cell1 vs Run1).
    for genome in sorted_genomes:
        prefix = f"{genome}_"
        mask = matched.isna() & seq_series.str.startswith(prefix)
        if mask.any():
            matched.loc[mask] = genome

    return matched


def load_external_label_sets(external_labels_path: str, og_id: str, seq_ids: list[str]) -> OrderedDict[str, pd.Series]:
    """Load genome- and proteome-level external labels aligned to sequence IDs.

    Accepted inputs for --external-labels:
    1) A directory containing:
       - achro_genome_labels.tsv (with a 'genome' column)
       - achro_proteome_labels.tsv (with 'og_id' and 'gene_id' columns)
    2) A single proteome labels TSV (legacy mode).
    """
    path = Path(external_labels_path)
    seq_index = pd.Index(seq_ids, name="sequence_id")
    out: OrderedDict[str, pd.Series] = OrderedDict()

    genome_df = None
    proteome_df = None

    if path.is_dir():
        genome_path = path / "achro_genome_labels.tsv"
        proteome_path = path / "achro_proteome_labels.tsv"
        if genome_path.exists():
            genome_df = pd.read_csv(genome_path, sep="\t", dtype="string")
        if proteome_path.exists():
            proteome_df = pd.read_csv(proteome_path, sep="\t", dtype="string")
        if genome_df is None and proteome_df is None:
            raise ValueError(
                "No label tables found in external-labels directory; expected achro_genome_labels.tsv and/or achro_proteome_labels.tsv"
            )
    else:
        single_df = pd.read_csv(path, sep="\t", dtype="string")
        if {"og_id", "gene_id"}.issubset(single_df.columns):
            proteome_df = single_df
        elif "genome" in single_df.columns:
            genome_df = single_df
        else:
            raise ValueError(
                "External labels TSV must be a proteome table (og_id,gene_id,...) or a genome table (genome,...)"
            )

    if genome_df is not None:
        if "genome" not in genome_df.columns:
            raise ValueError("Genome labels table must contain a 'genome' column")
        genome_df = genome_df.drop_duplicates(subset=["genome"], keep="first").set_index("genome")
        seq_genome = find_sequence_genome(seq_ids, genome_df.index)
        for col in genome_df.columns:
            values = seq_genome.map(genome_df[col])
            out[f"genome:{col}"] = values.fillna("NA").astype("string")

    if proteome_df is not None:
        if not {"og_id", "gene_id"}.issubset(proteome_df.columns):
            raise ValueError("Proteome labels table must contain 'og_id' and 'gene_id' columns")
        prot = proteome_df[proteome_df["og_id"] == og_id].copy()
        prot = prot.drop_duplicates(subset=["gene_id"], keep="first").set_index("gene_id")
        prot = prot.reindex(seq_index)

        for col in prot.columns:
            if col in {"og_id", "gene_id"}:
                continue
            out[f"proteome:{col}"] = prot[col].fillna("NA").astype("string")

    if not out:
        raise ValueError(f"No external labels available for OG {og_id}")

    return out


def choose_k(x: np.ndarray) -> tuple[int, float]:
    n = x.shape[0]
    if n < 3:
        return 1, float("nan")

    best_k = 2
    best_score = -1.0
    max_k = min(n - 1, 10)

    for k in range(2, max_k + 1):
        km = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = km.fit_predict(x)
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(x, labels)
        if score > best_score:
            best_score = score
            best_k = k

    return best_k, best_score


def main() -> int:
    args = parse_args()

    # Preserve the raw matrix (with NaNs) for plotting, but use a filled
    # copy for PCA and clustering computations.
    mat_raw = pd.read_csv(args.matrix, index_col=0)
    # Record NaN summary for this OG: number of sequences and rows with any NaN
    try:
        n_sequences = int(mat_raw.shape[0])
        rows_with_nan = int(mat_raw.isna().any(axis=1).sum())
        # Write summary into the repository's test_results/summary directory
        script_dir = Path(__file__).resolve().parent
        repo_root = script_dir.parent
        summary_path = repo_root / "test_results" / "summary" / "nan_rows_summary.csv"
        summary_dir.mkdir(parents=True, exist_ok=True)
        row = {
            "og_id": args.og_id,
            "matrix": args.matrix,
            "n_sequences": n_sequences,
            "rows_with_nan": rows_with_nan,
        }
        df_row = pd.DataFrame([row])
        write_header = not summary_path.exists()
        df_row.to_csv(summary_path, mode="a", header=write_header, index=False)
    except Exception:
        # Don't fail the clustering run for summary-writing issues.
        pass
    mat = mat_raw.copy()
    # Replace NaN distances with per-column means; fallback to global mean.
    mat = mat.apply(lambda c: c.fillna(c.mean()), axis=0)
    if mat.isna().values.any():
        global_mean = np.nanmean(mat.values)
        mat = mat.fillna(global_mean)

    x = mat.values.astype(float)
    seq_ids = mat.index.astype(str).tolist()

    n = x.shape[0]
    n_components = 2 if n >= 2 else 1
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(x)

    if n >= 3:
        k, sil = choose_k(x)
        if k == 1:
            labels = np.zeros(n, dtype=int)
        else:
            labels = KMeans(n_clusters=k, random_state=42, n_init=20).fit_predict(x)
    else:
        k, sil = 1, float("nan")
        labels = np.zeros(n, dtype=int)

    # Save per-sequence cluster assignments.
    cl_df = pd.DataFrame(
        {
            "og_id": args.og_id,
            "mode": args.mode,
            "sequence_id": seq_ids,
            "cluster": labels,
        }
    )
    cl_df.to_csv(args.out_clusters, index=False)

    # Permutation test for silhouette significance.
    sil_pval = silhouette_pvalue(x, labels, sil, n_perm=args.n_perm) if (k > 1 and not np.isnan(sil)) else float("nan")

    # Save model metadata.
    exp_var = pca.explained_variance_ratio_.tolist() if n_components > 1 else [1.0]
    meta_rows = [
        {
            "og_id": args.og_id,
            "mode": args.mode,
            "n_sequences": n,
            "n_clusters": int(k),
            "silhouette": sil,
            "silhouette_pvalue": sil_pval,
            "pc1_var": exp_var[0] if len(exp_var) > 0 else np.nan,
            "pc2_var": exp_var[1] if len(exp_var) > 1 else np.nan,
            "label_type": "internal",
            "label_name": "kmeans",
            "ari": float("nan"),
            "n_labels": len(set(labels)),
        }
    ]

    external_label_sets: OrderedDict[str, pd.Series] = OrderedDict()
    if args.external_labels is not None:
        external_label_sets = load_external_label_sets(args.external_labels, args.og_id, seq_ids)
        for full_label_name, ext_series in external_label_sets.items():
            label_type, label_name = full_label_name.split(":", 1)
            # Exclude missing/NA-coded values from ARI and silhouette calculations.
            # `load_external_label_sets` fills missing values with the string "NA",
            # so treat both pandas NA and the literal "NA" as missing here.
            valid_mask = ext_series.notna() & (ext_series != "NA")
            mask_arr = valid_mask.values
            valid_count = int(valid_mask.sum())

            if valid_count < 2:
                # Not enough valid observations to compute scores
                ari = float("nan")
                sil_ext = float("nan")
                sil_ext_pval = float("nan")
                n_ext = 0
            else:
                # Work with the subset that has non-missing labels
                ext_series_valid = ext_series[valid_mask]
                nums, names = map_categories_with_names(ext_series_valid)
                n_ext = len(names)

                if 2 <= n_ext < valid_count:
                    ari = adjusted_rand_score(labels[mask_arr], nums)
                    sil_ext = silhouette_score(x[mask_arr], nums)
                    sil_ext_pval = silhouette_pvalue(x[mask_arr], nums, sil_ext, n_perm=args.n_perm)
                else:
                    ari = float("nan")
                    sil_ext = float("nan")
                    sil_ext_pval = float("nan")

            meta_rows.append(
                {
                    "og_id": args.og_id,
                    "mode": args.mode,
                    "n_sequences": n,
                    "n_clusters": int(k),
                    "silhouette": sil_ext,
                    "silhouette_pvalue": sil_ext_pval,
                    "pc1_var": exp_var[0] if len(exp_var) > 0 else np.nan,
                    "pc2_var": exp_var[1] if len(exp_var) > 1 else np.nan,
                    "label_type": label_type,
                    "label_name": label_name,
                    "ari": ari,
                    "n_labels": int(n_ext),
                }
            )

    meta_df = pd.DataFrame(meta_rows)
    # Round numeric columns to 4 decimals for reproducible CSV outputs
    meta_df = meta_df.copy()
    meta_df = meta_df.round(4)
    meta_df.to_csv(args.out_meta, index=False)

    # Build multi-panel PCA plot.
    label_sets_num: dict[str, np.ndarray] = OrderedDict()
    label_names: dict[str, list[str]] = {}

    # Internal clustering: numeric labels; keep numeric names as strings
    label_sets_num["internal"] = labels
    label_names["internal"] = [str(i) for i in sorted(set(labels))]

    for full_label_name, ext_series in external_label_sets.items():
        nums, names = map_categories_with_names(ext_series)
        label_sets_num[full_label_name] = nums
        label_names[full_label_name] = names

    n_panels = len(label_sets_num)
    ncols = min(n_panels, 3)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    for ax, (full_label_name, label_values) in zip(axes.flat, label_sets_num.items()):
        # Use a discrete colormap sized to the number of unique labels
        names = label_names.get(full_label_name, [])
        n_labels = max(1, len(names))
        try:
            cmap = plt.cm.get_cmap("tab10", n_labels)
        except Exception:
            cmap = plt.cm.get_cmap("tab10")

        if n_components == 2:
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1], c=label_values, cmap=cmap, s=45, alpha=0.9, vmin=0, vmax=max(0, n_labels - 1)
            )
            ax.set_xlabel(f"PC1 ({exp_var[0] * 100:.1f}% var)")
            ax.set_ylabel(f"PC2 ({exp_var[1] * 100:.1f}% var)")
        else:
            scatter = ax.scatter(coords[:, 0], np.zeros_like(coords[:, 0]), c=label_values, cmap="tab10", s=45, alpha=0.9)
            ax.set_xlabel("PC1")
            ax.set_yticks([])

        if full_label_name == "internal":
            title = f"Internal clustering (k={k})"
        else:
            label_type, label_name = full_label_name.split(":", 1)
            title = f"External {label_type}: {label_name}"
        ax.set_title(title)
        cbar = fig.colorbar(scatter, ax=ax, label="Label", fraction=0.046, pad=0.04)
        # Show one tick per category and label them with the category names
        ticks = list(range(n_labels))
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(label_names.get(full_label_name, [str(t) for t in ticks]))
        cbar.ax.tick_params(size=0)

    for ax in axes.flat[n_panels:]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=150)
    plt.close(fig)

    # Heatmap of the distance matrix: use the raw matrix so original NaNs
    # are visible and will be masked (shown in grey).
    plot_heatmap(mat_raw, labels, args.og_id, args.mode, args.out_heatmap, args.out_sequence_order)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
