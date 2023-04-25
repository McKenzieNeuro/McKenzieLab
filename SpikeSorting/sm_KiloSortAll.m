% find files that need to be sorted
warning off
topDir = [ ...
    {'R:\STDP\STDP4'} ;...
    ];

ix= 1;
clear bb
for i = 1:length(topDir)
    ok = dir(topDir{i});
    ok = ok(cellfun(@any,{ok(1:end).isdir}));
    ok = ok(3:end);
    for j = 8:length(ok)
        
        cd(fullfile(ok(j).folder,ok(j).name))
        basename = bz_BasenameFromBasepath(ok(j).name);
        if exist([basename '.dat']) & exist([basename '.xml'])
            
            
            
            try
            savepath = sm_KiloSortWrapper('basepath',pwd,'config_version','McKenzie','intan_version','combined');
            end
        end
        
    end
    
end







