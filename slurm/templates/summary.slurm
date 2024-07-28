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
  -e {expt_name} \
  -f {fastq_dir} \
  -m {metadata} \
  -r {regions} \
  -p {pipeline} \
  --summary_only

