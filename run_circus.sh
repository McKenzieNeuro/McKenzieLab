module load miniconda3
module load matlab/R2022a


cd /carc/scratch/projects/mckenzie2016183/code/matlab
wait
source activate circus
wait


spyking-circus $1 -H $CARC_NODEFILE -c $SLURM_CPUS_PER_TASK 
wait
spyking-circus $1 -H $CARC_NODEFILE -m converting -c $SLURM_CPUS_PER_TASK 
#wait
#source deactivate
#wait
#matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "sm_reClusterKlustakwikCarc('$1'); exit;"


