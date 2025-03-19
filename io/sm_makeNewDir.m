function sm_makeNewDir(dirN)


if ~exist(fileparts(dirN))
    sm_makeNewDir(fileparts(dirN))
else
    mkdir(dirN)
    
end


if ~exist(dirN)
    sm_makeNewDir(dirN)
end
end