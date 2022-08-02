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
                [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd,'streams',{'x500D','x450D'},'isosbestic','x450D');
                
                streams{idxDA} = 'x500D';
            else
                %some were run on NE settings (COMPARE LATER)
                [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd);
                streams{idxDA} = 'x465A';
            end
            
            %smooth data
            signal_DFoF = nanconvn(signal_DFoF,k');
            
            %get indicded around reward
            ix = sm_getIndicesAroundEvent([data.epocs.U11_.onset;data.epocs.U12_.onset],10,10,fs,length(signal_DFoF));
            
            %save for DA
            DA{idxDA} = signal_DFoF(ix);
            subjID_DA{idxDA} = a(isDA(1):isDA(1)+4);
            date_DA{idxDA} = data.info.date;
            
            powerDA(idxDA,:) = data.scalars.Fi1i.data([3 7],1);
            idxDA = idxDA+1;
        elseif ~isempty(isNE)
            %load iso corrected for NE (default settings)
            [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd);
            
            %smooth
            signal_DFoF = nanconvn(signal_DFoF,k');
            
            %get indices around reward
            ix = sm_getIndicesAroundEvent([data.epocs.U11_.onset;data.epocs.U12_.onset],10,10,fs,length(signal_DFoF));
            
            %save for NE
            NE{idxNE} = signal_DFoF(ix);
            date_NE{idxNE} = data.info.date;
            subjID_NE{idxNE} = a(isNE(1):isNE(1)+4);
            powerNE(idxNE,:) = data.scalars.Fi1i.data([3 7],1);
            idxNE = idxNE+1;
        end
        
    end
end