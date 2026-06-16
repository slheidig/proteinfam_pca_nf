#!/usr/bin/env python3
"""Split merged b2b TSV back into per-OG TSV files by OG_ID|sequence_id.

Some b2bTools outputs normalize the separator to an underscore, so we also
accept OG_ID_sequence_id and recover the OG ID from the first separator.
"""

from __future__ import annotations

import argparse

try:
    import pandas as pd
except ImportError:  # pragma: no cover - container should normally provide pandas
    pd = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("merged_tsv")
    p.add_argument("names_map", nargs="?", default=None,
                   help="optional TSV with columns: OG\toriginal_name\tstripped_name")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if pd is None:
        raise ImportError("pandas is required for reattach_b2b.py")

    df = pd.read_csv(args.merged_tsv, sep="\t")
    if df.empty:
        return 0
    if "sequence_id" not in df.columns:
        raise ValueError("Input TSV missing required column: sequence_id")

    extracted = df["sequence_id"].astype("string").str.extract(
        r"^(OG[^|_]+)[|_](.+)$"
    )
    extracted.columns = ["og_id", "original_id"]
    valid_rows = extracted["og_id"].notna() & extracted["original_id"].notna()
    df = df.loc[valid_rows].copy()
    if df.empty:
        return 0

    df["sequence_id"] = extracted.loc[valid_rows, "original_id"].to_numpy()
    df.insert(0, "og_id", extracted.loc[valid_rows, "og_id"].to_numpy())

    # If a names mapping is provided, restore original sequence names
    # where b2bTools removed characters (e.g., dots).
    if args.names_map:
        map_df = pd.read_csv(args.names_map, sep="\t", dtype=str)
        if not set(["OG", "original_name", "stripped_name"]).issubset(map_df.columns):
            raise ValueError("names_map missing required columns: OG, original_name, stripped_name")
        map_df = map_df.rename(columns={"OG": "og_id", "original_name": "orig_name", "stripped_name": "stripped_name"})
        # merge on og_id and the stripped name (current sequence_id)
        df = df.merge(
            map_df[["og_id", "stripped_name", "orig_name"]],
            how="left",
            left_on=["og_id", "sequence_id"],
            right_on=["og_id", "stripped_name"],
        )
        # where mapping exists, replace sequence_id with the original name
        df["sequence_id"] = df["orig_name"].where(df["orig_name"].notna(), df["sequence_id"])
        df = df.drop(columns=[c for c in ("stripped_name", "orig_name") if c in df.columns])

    for og_id, group in df.groupby("og_id", sort=False):
        group.drop(columns=["og_id"]).to_csv(
            f"{og_id}_b2b.tsv",
            sep="\t",
            index=False,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
