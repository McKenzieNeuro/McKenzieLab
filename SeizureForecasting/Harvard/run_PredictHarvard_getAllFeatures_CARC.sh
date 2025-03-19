#!/bin/bash
#SBATCH --time=7-00:00:00       # walltime
#SBATCH --partition=neuro-hsc
#SBATCH -J "EMD"       # job name
#SBATCH --array=1-999# job array of input size
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks-per-node=1

module load matlab/R2021a



cd /carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/code/matlab

matlab -singleCompThread -nodisplay  -nodesktop  -r "sm_PredictHarvard_getAllFeatures_CARC($SLURM_ARRAY_TASK_ID); exit;"


