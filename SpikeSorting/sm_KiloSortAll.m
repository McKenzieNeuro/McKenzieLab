% find files that need to be sorted
warning off
topDir = [ ...
    {'E:\data\Paudel'} ;...
    ];

ix= 1;
clear bb
for i = 1:length(topDir)
    ok = dir(topDir{i});
    ok = ok(cellfun(@any,{ok(1:end).isdir}));
    ok = ok(3:end);
    for j = 4:length(ok)
        
        subdir = [topDir{1} filesep ok(j).name];
        
        ok1 = dir(subdir);
        ok1 = ok1(cellfun(@any,{ok1(1:end).isdir}));
        ok1 = ok1(3:end);
        kp = ~ (cellfun(@any,regexpi({ok1.name}','histology')) | cellfun(@any,regexpi({ok1.name}','impedance')));
        ok1 = ok1(kp);
        %exclude histology and impedance
        
        
        for k = 4:length(ok1)
            
            subsubdir = [subdir filesep ok1(k).name];
            
            ok2 = dir(subsubdir);
            
            %check if kilosort ran
            disp(['working on ' subsubdir])
            if ~any(cellfun(@any,regexp({ok2.name},'kilosortDone')))
                
                try
                    dirName  = subsubdir;
                    
                    basename = bz_BasenameFromBasepath(dirName);
                    
                 %   sm_Process_Intan_kilosort(basename,dirName,'intan_version','Combined');
                    
                    
                    
                    
                    
              savepath = sm_KiloSortWrapper('basepath',dirName,'config_version','McKenzie','intan_version','combined');
                    
                end
                
            end
            
            
            
            
        end
        
        
        
    end
end