#!/bin/bash

analyse_jid=$(sbatch --parsable {run_dir}/analyse.slurm)
sbatch --parsable --dependency=afterany:$analyse_jid {run_dir}/summary.slurm
