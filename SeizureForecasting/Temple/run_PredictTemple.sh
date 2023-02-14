#!/bin/bash
###SBATCH --time=010:00:00       # walltime
#SBATCH --partition=neuro-hsc
#SBATCH -J "PredictTemple1"       # job name
#SBATCH --array=1          # job array of input size
#SBATCH --nodes=1
#SBATCH --tasks=1


module load matlab/R2022a



cd /carc/scratch/projects/mckenzie2016183/code/matlab

matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "sm_PredictTemple_getAllFeatures_CARC($SLURM_ARRAY_TASK_ID); exit;"