fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);


kp = cellfun(@any,regexp(fils,'Novel Env')) & cellfun(@any,regexpi(fils,'transition'));
fils = fils(kp);

[dirs,b] = fileparts(fils);
dirs = unique(dirs);
%%
tot = 0;
i=3;
   cd(dirs{i})
   
   
 
    
%%




  sm_getArenaEdges(dirs{i})

  
  [dirs,b] = fileparts(fils);