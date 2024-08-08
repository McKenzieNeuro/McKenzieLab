% find files that need to be sorted
warning off
% topDir = [ ...
%     {'R:\STDP\STDP6'} ;...
%     {'R:\STDP\STDP8'} ;...
%     {'R:\STDP\STDP9'} ;...
%     {'R:\STDP\STDP12'} ;...
%     ];

topDir{1} = 'R:\WSun\Bicon Behav\Ctxdis2';
ix= 1;
clear bb
for i = 1:length(topDir)
    ok = dir(topDir{i});
    ok = ok(cellfun(@any,{ok(1:end).isdir}));
    ok = ok(3:end);
    for j = 16:length(ok)
        
        cd(fullfile(ok(j).folder,ok(j).name))
        basename = 'amplifier_analogin_auxiliary_int16';
        fils = getAllExtFiles(pwd,'dat',1);
        fils = fils(contains(fils,'amplifier_analogin_auxiliary_int16.dat'));
        
         fils1 = getAllExtFiles(pwd,'mat',0);
        fils1 = fils1(contains(fils1,'kilosortDone'));
        
        
        if length(fils)==1 && isempty(fils1)
            
            
            
            try
            dirN = fileparts(fils{1});
                savepath = sm_KiloSortWrapper('basepath',dirN,'basename',basename,'config_version','McKenzie','intan_version','combined');
            end
        end
        
    end
    
end







