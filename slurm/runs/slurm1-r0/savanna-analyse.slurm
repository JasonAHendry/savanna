#!/bin/bash -l
#SBATCH --job-name=analyse
#SBATCH --output=logs/analyse-jid%A-%a.out
#SBATCH --error=logs/analyse-jid%A-%a.out
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array 9,10,12,13
#SBATCH --mem=16GB
#SBATCH --time=23:00:00

conda activate savanna

savanna analyse \
  -e slurm1 \
  -m <savanna.util.metadata.MetadataTableParser object at 0x7f7c7820bf90> \
  -r example_data/ag_set1/metadata/nomadsSID.amplicons.bed \
  -p ento \
  --barcode $SLURM_ARRAY_TASK_ID

