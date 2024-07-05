#!/bin/bash

basecall_jid=$(sbatch --parsable {run_dir}/basecall.slurm)
analyse_jid=$(sbatch --parsable --dependency=afterok:$basecall_jid {run_dir}/analyse.slurm)
sbatch --parsable --dependency=afterany:$analyse_jid {run_dir}/summary.slurm
