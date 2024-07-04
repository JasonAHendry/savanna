#!/bin/bash -l
#SBATCH --job-name=basecall
#SBATCH --output=logs/basecall-jid%A-%a.out
#SBATCH --error=logs/basecall-jid%A-%a.err
#SBATCH --chdir=./
#SBATCH --ntasks=1
# --------------------
# GPU Settings
# --------------------
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#SBATCH --time=24:00:00

conda activate savanna

savanna basecall \
 -e {expt_name} \
 -p {pod5_dir} \
 -m {model} \
 -q {min_qscore}

savanna demultiplex \
 -e {expt_name} \
 -k {kit} \
 -t guppy \
 --strict
