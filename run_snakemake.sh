#!/bin/bash -l
#SBATCH --job-name=snakemake_grn
#SBATCH --account=upzenk
#SBATCH --partition=standard
#SBATCH --qos=serial
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --output=/home/jiafzhou/GRN_org/logs/snakemake_driver_%j.out
#SBATCH --error=/home/jiafzhou/GRN_org/logs/snakemake_driver_%j.err

set -euo pipefail

source /work/upzenk/miniconda3/bin/activate
cd /home/jiafzhou/GRN_org

conda activate snakemake

snakemake \
  -s workflow/Snakefile \
  --profile profiles/slurm \
  --jobs 5 \
  --restart-times 0 \
  all