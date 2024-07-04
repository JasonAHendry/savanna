#!/bin/bash

basecall_jid=$(sbatch --parsable {run_dir}/savanna-basecall.slurm)
analyse_jid=$(sbatch --parsable --dependency=afterok:$basecall_jid {run_dir}/savanna-analyse.slurm)
sbatch --parsable --dependency=afterany:$analyse_jid {run_dir}/savanna-summary.slurm
