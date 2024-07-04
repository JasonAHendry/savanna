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
 -e slurm1 \
 -p None \
 -m sup \
 -q 10

savanna demultiplex \
 -e slurm1 \
 -k SQK-NBD114-96 \
 -t guppy \
 --strict
