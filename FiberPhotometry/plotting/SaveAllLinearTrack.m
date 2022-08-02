fils = getAllExtFiles(pwd,'mat',1);

%find all linear track files
kp = cellfun(@any,regexpi(fils,'trackON'));

fils = fils(kp);

%%
k  = gaussian2Dfilter([10000 1],100);
idxDA = 1;
idxNE = 1;
clear DA NE powerDA powerNE date_NE date_DA subjID_NE subjID_DA

%loop over all files
for i = 1:length(fils)
    
    %load in fiber photometry data
    [a,b]  =fileparts(fils{i});
    cd(a)
    data = TDTbin2mat(a);
    
    %decide if it is a NE or DA animal (note DACSD1 has a different naming
    %convention
    
    isNE = regexp(a,'NE2');
    isDA = regexp(a,'DA2');
    
    if isempty(isDA) & isempty(isNE)
        isDA = regexp(a,'DACSD');
    end
    
    %detect if the sync pulses were set up
    if isfield(data.epocs,'U11_') && isfield(data.epocs,'U12_')
        
        if ~isempty(isDA)
            
            %for the DA animals, check if the right LEDs were used
            if isfield(data.streams,'x500D')
                [GEFI.corrected ,GEFI.ts_data,GEFI.fs] = sm_getSignal_DFoF(pwd,'streams',{'x500D','x450D'},'isosbestic','x450D');
                GEFI.signal = 'x500D';
                GEFI.isosbestic = 'x450D';
                
            else
                %some were run on NE settings (COMPARE LATER)
                [GEFI.corrected,GEFI.ts_data,GEFI.fs] = sm_getSignal_DFoF(pwd);
                GEFI.signal = 'x465A';
                GEFI.isosbestic = 'x405A';
            end
            
            
            
            
            
            
            GEFI.subjName =  a(isDA(1):isDA(1)+4);
            
            
            GEFI.indicator = 'DA';
            
        elseif ~isempty(isNE)
            %load iso corrected for NE (default settings)
            [GEFI.corrected,GEFI.ts_data,GEFI.fs] = sm_getSignal_DFoF(pwd);
            
            
            GEFI.signal = 'x465A';
            GEFI.isosbestic = 'x405A';
            GEFI.subjName = a(isNE(1):isNE(1)+4);
            
            
        end
        
        GEFI.indicator = 'NE';
        GEFI.reward1 = data.epocs.U11_.onset;
        GEFI.reward2 = data.epocs.U12_.onset;
        GEFI.date = data.info.date;
        GEFI.power_signal  = data.scalars.Fi1i.data([3],1); %check this
        GEFI.power_isos  = data.scalars.Fi1i.data([7],1);%check this
        dash = regexp(fils{i},'-');
        fdate = fils{i}(dash(1)+1:dash(2)-1);
        fname = [GEFI.subjName '_' fdate '.GEFI.mat'];
        save(fname,'GEFI')
        
        clear GEFI
    end
end