fils = getAllExtFiles('R:\DANEHippocampalResponse','csv',1);
fils = fils(cellfun(@any,regexpi(fils,'linear')));

[dirs,b] = fileparts(fils);
dirs = unique(dirs);
%%
for i = 1:length(dirs)
    try
        disp(['working on ' dirs{i}])
        saveAllNeuralPosition(dirs{i})
        
    catch
        disp('bad dir')
        
    end
    i
end



