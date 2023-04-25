#!/bin/bash
#SBATCH --partition=neuro-hsc
#SBATCH -J "spyking-circus"       # job name
#SBATCH --array=1           # job array of input size
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8


datadir="/carc/scratch/projects/mckenzie2016183/data/spikeSorting/STDP/STDP/STDP4"


module load parallel


find $datadir -name "*.dat" | parallel --jobs $SLURM_NTASKS --joblog $SLURM_JOB_NAME.joblog --resume /carc/scratch/projects/mckenzie2016183/code/BASH/run_circus.sh {}
#find $datadir -name "*.dat" | cat