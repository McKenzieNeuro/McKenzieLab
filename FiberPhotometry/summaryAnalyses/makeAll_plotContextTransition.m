fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);


%%
clear labels
for i = 1:length(fils)
    
   v = load(fils{i}); 
   labels{i}=v.data(:,1);
end
%%

fils(cellfun(@(a) any(cellfun(@any,regexp(a,'context2'))),labels))
