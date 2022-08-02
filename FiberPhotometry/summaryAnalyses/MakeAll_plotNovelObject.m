

%%
fils = getAllExtFiles('R:\DANEHippocampalResponse','Tbk',1);
kp = (cellfun(@any,regexpi(fils,'SOR')));
fils = fils(kp);
[directory] = fileparts(fils);
%%
nSessions = length(directory);
%%
for i = 1:nSessions
cd(directory{i})


end


for i = 1:nSessions


[left{i},right{i},ts_PETH{i}]  = plotNovelObject(directory{i});
i

end