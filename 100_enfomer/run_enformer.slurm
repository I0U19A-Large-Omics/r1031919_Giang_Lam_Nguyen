#!/bin/bash

#SBATCH --job-name=fake_enformer

#SBATCH --account=lp_edu_large_omics

#SBATCH --clusters=wice

#SBATCH --time=00:05:00

#SBATCH --cpus-per-task=1

#SBATCH --mem=1G

#SBATCH --output=enformer.%j.out

#SBATCH --error=enformer.%j.err



set -euo pipefail



cd /data/leuven/386/vsc38681/r1031919_Giang_Lam_Nguyen/100_enfomer

export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2026/bin:$PATH



python score_enformer.py \

  ../030_snakemake_snpcalling/50.annotated/annotated_snps.vcf \

  --output enformer_scores.tsv \

  --threads "${SLURM_CPUS_PER_TASK}"



head -100 enformer_scores.tsv > enformer_scores.head100.tsv

