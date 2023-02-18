#!/bin/bash
#SBATCH --time=10:00:00       # walltime
#SBATCH --partition=neuro-hsc
#SBATCH -J "spyking-circus"       # job name
#SBATCH --array=1           # job array of input size
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8


datadir="/carc/scratch/projects/mckenzie2016183/data/spikeSorting/spikeDemo"


module load parallel
module load miniconda3
module load matlab/R2022a

find $datadir -name "*.dat" | parallel --jobs $SLURM_NTASKS --joblog $SLURM_JOB_NAME.joblog --resume /carc/scratch/projects/mckenzie2016183/code/BASH/run_circus.sh {}