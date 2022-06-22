#!/bin/bash -login
#SBATCH -t 1:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1gb

set -o nounset
set -o errexit
set -x

echo hello, world

# make a directory specific to user and job
export MYTMP=/scratch/${USER}/slurm_${SLURM_JOBID}
mkdir -p $MYTMP

# force clean it up
function cleanup() { rm -rf $MYTMP; }
trap cleanup EXIT

# <run your code, telling it to use that dir;
echo running in $(pwd)

echo to copy results out, you could do
echo cp $MYTMP/files_you_want_to_keep $SLURM_SUBMIT_DIR

echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
