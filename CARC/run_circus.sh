#!/bin/bash



file=$1

cd /carc/scratch/projects/mckenzie2016183/code

source activate circus



spyking-circus $file -H $CARC_NODEFILE -c $SLURM_CPUS_PER_TASK 
spyking-circus $file -H $CARC_NODEFILE -m converting -c $SLURM_CPUS_PER_TASK 
matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "sm_reClusterKlustakwikCarc('$file'); exit;"


