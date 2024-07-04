#!/bin/bash

basecall_jid=$(sbatch --parsable slurm/runs/slurm1-r3/savanna-basecall.slurm)
analyse_jid=$(sbatch --parsable --dependency=afterok:$basecall_jid slurm/runs/slurm1-r3/savanna-analyse.slurm)
sbatch --parsable --dependency=afterany:$analyse_jid slurm/runs/slurm1-r3/savanna-summary.slurm
