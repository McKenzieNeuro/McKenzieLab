#!/bin/bash
###SBATCH --time=010:00:00       # walltime
#SBATCH --partition=neuro-hsc
#SBATCH -J "addOne"       # job name
#SBATCH --array=1-1000          # job array of input size
#SBATCH --nodes=1
#SBATCH --tasks=1


file=$1

module load matlab/R2022a



cd /carc/scratch/projects/mckenzie2016183/code

matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "addOne('$file',$SLURM_ARRAY_TASK_ID); exit;"


