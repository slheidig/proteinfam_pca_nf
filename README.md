# og-b2bpca: How-To Guide

This pipeline detects internal structural-dynamic clusters in orthogroups using:
- b2bTools backbone predictions
- MAFFT multiple sequence alignment
- MMseqs2 pairwise alignment
- Distance matrices + PCA + silhouette-based KMeans clustering

## 1) What You Need

- Nextflow (>= 24.04)
- One container profile:
  - `-profile singularity` (Apptainer/Singularity)
  - or `-profile docker`
- Input directory with orthogroup FASTA files named `*.fa`

Each FASTA file is treated as one orthogroup.

## 2) Quick Start (test dataset)

From the pipeline root:

```bash
module load Nextflow
nextflow run . -profile test,singularity -c tests/nextflow.config
```

This runs a small built-in test and writes outputs to `test_results/`.

## 3) Run on Your Own Data

```bash
module load Nextflow
nextflow run . \
  -profile singularity \
  --og_dir /path/to/Orthogroup_Sequences \
  --outdir results \
  --min_seqs 10 \
  --b2b_batch_size 1000
```

Docker alternative:

```bash
nextflow run . \
  -profile docker \
  --og_dir /path/to/Orthogroup_Sequences \
  --outdir results
```

## 4) Key Parameters

- `--og_dir` (required): folder containing orthogroup FASTA files (`*.fa`)
- `--outdir` (default `results`): output root folder
- `--min_seqs` (default `10`): skip OGs below this sequence count
- `--b2b_batch_size` (default `1000`): batch size for b2bTools prediction stage

Show CLI help:

```bash
nextflow run . --help
```

## 5) What the Pipeline Produces

Inside `--outdir`:

- `b2b_predictions/`
  - per-OG `*_b2b.tsv`
- `mafft/<OG_ID>/`
  - MAFFT alignment `*.aln.fa`
- `mmseqs2/<OG_ID>/`
  - MMseqs2 pairwise alignment TSV
- `distance_matrices/<OG_ID>/`
  - `*_b2b_dist.csv` distance matrix
  - `*_pair_meta.tsv` per-pair comparison metadata
- `pca/<OG_ID>/`
  - `*_pca.png` PCA scatter colored by cluster
  - `*_heatmap.png` distance heatmap with sequence labels
  - `*_clusters.csv` sequence-to-cluster assignments
  - `*_pca_meta.csv` clustering metadata (including silhouette and permutation p-value)
- `summary/`
  - `cluster_counts_summary.csv` per-OG internal cluster count
  - `cluster_count_histogram.csv` binned histogram table (`1..10` and `10+`)
  - `cluster_counts_barplot.png` histogram plot of internal cluster counts per mode

## 6) Interpreting the Main Outputs

- `*_clusters.csv`: final cluster label for each sequence.
- `*_pca_meta.csv`:
  - `n_clusters`: selected KMeans cluster count
  - `silhouette`: separation quality (higher is better)
  - `silhouette_pvalue`: permutation significance estimate
- `*_heatmap.png`: matrix sorted by cluster; red boundaries indicate cluster transitions.
- `cluster_count_histogram.csv`: number of OGs per internal-cluster bin and mode.

## 7) Useful Rerun Patterns

Resume from cache after interruptions:

```bash
nextflow run . -profile singularity --og_dir /path/to/Orthogroup_Sequences --outdir results -resume
```

Use a custom config on HPC:

```bash
nextflow run . -profile slurm,singularity -c your_config.config --og_dir /path/to/Orthogroup_Sequences --outdir results
```
