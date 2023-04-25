
% collect all data
topdir = 'R:\DGregg\NeuralData\LTD 10.0';
files1 = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files1,'sessiondata');

files1 = files1(keepFiles);

topdir = 'R:\DGregg\NeuralData\EDS';
files2 = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files2,'sessiondata');

files2 = files2(keepFiles);

files = [files1;files2];
chReg = [...
    {'PrL(L)' };...
    {'PrL(R)' };...
    {'AVT(L)' };...
    {'BLA(R)' };...
    {'CA1(L)' };...
    {'CA1(R)' };...
    {'gRSC (L)'};...
    {'LDT(L1)'}];
%%

%u_lfp_ied_all = nan(2500,8,8,length(files));
for i = 97:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    if isfield(sessiondata,'IED')
        if ~iscell(sessiondata.lfpFile)
            % get all IEDs <stim1
            if ~isempty(sessiondata.stimON) && ~isnan(sessiondata.stimON(1))
                IED = cellfun(@(a) a(a<sessiondata.stimON(1)),sessiondata.IED,'uni',0);
            else
                IED = sessiondata.IED;
            end
            xml = LoadXml([sessiondata.lfpFile(1:end-3) 'xml']);
            fileDur = sessiondata.ts(end);
            lfpfil = sessiondata.lfpFile;
        else
            
            IED = cellfun(@(a) a(a<  sessiondata.fileDur(1)),sessiondata.IED,'uni',0);
            xml = LoadXml([sessiondata.lfpFile{1}(1:end-3) 'xml']);
            lfpfil = sessiondata.lfpFile{1};
            fileDur = sessiondata.fileDur(1);
        end
        [kp_ch,chIDX] = ismember(sessiondata.channel,chReg);
        % loop through each IED
        
        for j = 1:length(IED)
            
            if chIDX(j)>0
                all_dat = nan(2500,8,length(IED{j}));
                for k = 1:length(IED{j})
                    
                    ts = IED{j}(k);
                    
                    if ts>1 & ts <fileDur-1
                        d = LoadBinary(lfpfil,'nchannels',xml.nChannels,'frequency',xml.lfpSampleRate,'channels',sessiondata.channelID,...
                            'start',ts-1,'duration',2);
                        all_dat(:,chIDX(chIDX>0),k) = d(:,kp_ch);
                    end
                end
                
                  u_lfp_ied(:,:,chIDX(j)) = nanmean(all_dat,3);
            end
            
          
        end
        
        u_lfp_ied_all(:,:,:,i) = u_lfp_ied;
        i
    end
end

%%
close all
for i = [1 3 5 8]
    
figure

for j =  [1:6 8]
plot((-1250:1250-1)/1250,.195*(nanmedian(u_lfp_ied_all(:,j,i,:),4)-j*1000),'k')
hold on
end
title(chReg{i})
xlim([-.3 .3])
end