fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);


kp = cellfun(@any,regexp(fils,'Novel Env')) & cellfun(@any,regexpi(fils,'transition'));
fils = fils(kp);

[dirs,b] = fileparts(fils);
dirs = unique(dirs);
%%
tot = 0;
for i = 1:length(dirs)
   cd(dirs{i})
   
   
    if ~exist('arena_edges.mat')
    
        
    
     error('work here')
     
        
        
    end
   i
end
   tot 

   
%%




  sm_getArenaEdges(dirs{i})

  
  [dirs,b] = fileparts(fils);