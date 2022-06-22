#!/bin/bash -login
#SBATCH -t 1:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"
conda activate sgc

set -o nounset
set -o errexit
set -x

echo hello, world

echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
