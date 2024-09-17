function sm_ConcatDats(fnamesin,fnameout)
%this function merges fnamesin and copies to fnameout
%fnamesin is a cell array with full file paths
cmd = 'COPY /B';

for i = 1:length(fnamesin)
    
    if i ==1
        cmd = [cmd ' ' fnamesin{i} ];
    else
    cmd = [cmd ' + ' fnamesin{i} ];
    end
end

cmd = [cmd ' ' fnameout];

system(cmd)

end

