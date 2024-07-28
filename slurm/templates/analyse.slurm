#!/bin/bash -l
#SBATCH --job-name=analyse
#SBATCH --output=logs/analyse-jid%A-%a.out
#SBATCH --error=logs/analyse-jid%A-%a.out
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array {array_str}
#SBATCH --mem=16GB
#SBATCH --time=23:00:00

conda activate savanna

savanna analyse \
  -e {expt_name} \
  -f {fastq_dir} \
  -m {metadata} \
  -r {regions} \
  -p {pipeline} \
  --barcode $SLURM_ARRAY_TASK_ID

