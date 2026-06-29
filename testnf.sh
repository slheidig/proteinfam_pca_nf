#!/bin/bash
#SBATCH --job-name=nf-og_b2b_pca
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G


APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer
export APPTAINER_CACHEDIR=$apptainercache

module load Nextflow

nextflow run $VSC_SCRATCH_VO_USER/achromatium/35_b2b_pca/og-b2bpca/main.nf \
    -profile slurm,singularity -resume \
    --external_labels $VSC_SCRATCH_VO_USER/achromatium/35_b2b_pca/labels \
    --og_dir $VSC_SCRATCH_VO_USER/achromatium/35_b2b_pca/og-b2bpca/tests/data/ogs \
    --apptainer_cache $APPTAINERCACHE \
    --distance_type seq_distance \
    --outdir $VSC_SCRATCH_VO_USER/achromatium/35_b2b_pca/tests/seqmatrix 
    