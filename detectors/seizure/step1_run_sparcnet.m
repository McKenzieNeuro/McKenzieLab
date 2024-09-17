clc;close all;
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call python model sparcnet in env "sparcnet" (see config instructions how-to)
% alternatively run following w/o matlab
% $ conda activate sparcnet
% $ python ./callbacks/sparcnet1.0/fcn_run_sparcnet1.py *input_dir *output_dir *sampling_rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py_exe='/home/exx/anaconda3/envs/sparcnet/bin/python'; % $ which python
py_code='./callbacks/sparcnet1.0/fcn_run_sparcnet1.py';
    
input_dir='./data/raw/';
output_dir='./data/sparcnet/';
sampling_rate=200;

cmd=[py_exe,' ',[py_code,' ',input_dir,' ',output_dir,' ',num2str(sampling_rate)]];
system(cmd);