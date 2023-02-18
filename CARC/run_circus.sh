cd /carc/scratch/projects/mckenzie2016183/code/matlab
source activate circus



#spyking-circus $1 -H $CARC_NODEFILE -c $SLURM_CPUS_PER_TASK 
#spyking-circus $1 -H $CARC_NODEFILE -m converting -c $SLURM_CPUS_PER_TASK 

#source deactivate

matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "sm_reClusterKlustakwikCarc('$1'); exit;"


