#!/bin/bash
#SBATCH --time=10:00:00       # walltime
#SBATCH --partition=neuro-hsc
#SBATCH -J "spyking-circus"       # job name
#SBATCH --array=1           # job array of input size
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8


datadir="/users/mckenzie/data/spikeSorting"


module load parallel
module load miniconda3
module load matlab/R2022a

find $datadir -name "*int16.dat" | parallel --jobs $SLURM_NTASKS --joblog $SLURM_JOB_NAME.joblog --resume spyking-circus {} -H $CARC_NODEFILE -c $SLURM_CPUS_PER_TASK 