
module load matlab/R2022a


cd /carc/scratch/projects/mckenzie2016183/code/matlab
wait
matlab -singleCompThread -nodisplay  -nodesktop  -nojvm -r "sm_reClusterKlustakwikCarc('$1'); exit;"


