#!/bin/bash -l
#SBATCH --job-name=summary
#SBATCH --output=logs/summary-jid%A-%a.out
#SBATCH --error=logs/summary-jid%A-%a.out
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=23:00:00

conda activate savanna

savanna analyse \
  -e slurm1 \
  -m example_data/ag_set1/metadata/sample_info.csv \
  -r example_data/ag_set1/metadata/nomadsSID.amplicons.bed \
  -p plasmo \
  --summary_only

