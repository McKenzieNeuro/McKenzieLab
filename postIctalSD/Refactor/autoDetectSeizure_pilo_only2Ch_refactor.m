 % topDirs{1} = 'R:\DGregg\NeuralData\PTP';
 % topDirs{2} = 'R:\DGregg\NeuralData\PCP';
%topDirs{1} = 'R:\ASommer\PilocarpineRecordings\PTP_5.1';
topDirs{1} = 'R:\ASommer\PilocarpineRecordings\Tests\UnpluggedTest\';
% topDirs{1} = 'R:\DGregg\NeuralData\PTP';
fils = [];
for i = 1:length(topDirs)
    
    tmp = getAllExtFiles(topDirs{i},'mat',1);
    
    kp = contains(tmp,'sessiondata');
    fils = [fils;tmp(kp)];
end

%%
% 
% fils = [];
% for i = 1:length(topDirs)
% 
%     tmp = getAllExtFiles(topDirs{i},'szr',1);
% 
%     kp = contains(tmp,'evt');
%     fils = [fils;tmp(kp)];
% end

%%
% all_sz = [];alldiff=[];
% %get all seizures
% subj = '4.0';
% for i = 1:length(fils)
% 
%     dirN = fileparts(fils{i});
%     cd(dirN);
%     ev = LoadEvents(fils{i});
% 
%     if ~isempty(ev.time)
%         kp = (contains(lower(ev.description),'sz_on') |  contains(lower(ev.description),'seizurestart'))&contains(lower(ev.description),subj);
%         sz_on  = ev.time(kp);
%         kp = (contains(lower(ev.description),'sz_off') |  contains(lower(ev.description),'seizurestop'))&contains(lower(ev.description),subj);
% 
%         sz_off  = ev.time(kp);
% 
%         sz_ep = [sz_on sz_off];
% 
%         if ~isempty(sz_ep)
%             alldiff = [alldiff;diff(sz_ep(:,1))];
%             %find sessiondata
%             fil = getAllExtFiles(pwd,'mat',0);
% 
%             kp = contains(fil,'sessiondata') & contains(fil,subj);
% 
%             fil = fil(kp);
% 
%             if ~isempty(fil)
%                 load(fil{1})
%             else
%                 error('make session info file')
%             end
% 
% 
%             if any(kp)
%                 xml = sessiondata.xml;
% 
%                 for j = 1:size(sz_ep,1)
% 
%                     dur = sz_ep(j,2)-sz_ep(j,1);
% 
%                     d = LoadBinary(sessiondata.lfp_file,'nchannels',xml.nChannels,'channels', ...
%                         sessiondata.channelID,'frequency',sessiondata.xml.lfpSampleRate, ...
%                         'start',sz_ep(j,1)-60,'duration',120);
% 
%                     all_sz = cat(3,all_sz,d);
%                 end
%             end
% 
%         end
%     end
% 
% end
% 
% %%
% if ~isempty(all_sz)
% for i = 1:size(all_sz,3)
% 
%    figure 
%     plot(all_sz(:,5,i))
%     waitforbuttonpress
%     close all
% 
% 
% end
% end

%%

subjID =  'Test';

kp = contains(fils,subjID);

fils1 = fils(kp);

%%


FeatureFileOutput = 'R:\Analysis\McKenzieLab\postIctalSD\Refactor\Features\features_pilo18FEB26.mat';
load(FeatureFileOutput)
ClassifierFileOutput = 'R:\Analysis\McKenzieLab\postIctalSD\Refactor\Models\model_pilo_bilateral_18FEB26.mat';
load(ClassifierFileOutput)
%%

% make LFP file

for i = 1:length(fils1)
    clear sessiondata
    dirN = fileparts(fils1{i});
    fnamei = [dirN filesep 'amplifier.dat'];
    fnameo = [dirN filesep 'amplifier.lfp'];
    fxml = strrep(fnameo,'lfp','xml');
    try
        if ~exist(fnameo)
            bz_LFPfromDat(dirN,'basename','amplifier')
        end
    end
    load(fils1{i})
    sessiondata.lfp_file = fnameo;
    
    
    xml = LoadXml(fxml);
    sessiondata.xml = xml;
    
    save(fils1{i},'sessiondata');
    i
end
%%
% detect all subjID seizures


for i = 1:length(fils1)
    fout  = [fileparts(fils1{i}) filesep 'auto_sz_classifier_refactor.mat'];
   % if ~exist(fout)
        v = load(fils1{i});
        fileDur = sm_getFileDur(v.sessiondata.lfp_file);
        %get prediction
        estimateLabel =[];
        dat1 =[];
        % REMOVE BELOW THIS LATER 
            v.sessiondata.channelID = [1,5];
        % REMOVE ABOVE THIS LATER
        for tim = 1:floor(fileDur)-ops.durFeat
            data = LoadBinary(v.sessiondata.lfp_file,'nchannels',v.sessiondata.xml.nChannels,'frequency', ...
                v.sessiondata.xml.lfpSampleRate,'channels',v.sessiondata.channelID,'duration', ops.durFeat,'start',tim);
            
            
            features = sm_GetDataFeature_rat(data,tim,ops);
            
            dat1 = [dat1;features];
            
            if mod(tim,100)== 0
                [outpred,conf] = predict(rusTree,dat1);
                estimateLabel = [estimateLabel;outpred conf];
                dat1 =[];
            elseif tim > (fileDur- mod(fileDur,100))
                [outpred,conf] = predict(rusTree,dat1);
                estimateLabel = [estimateLabel;outpred conf];
                dat1 =[];
            end
            
        end
        sz_ep  =[find(diff([0;estimateLabel(:,1)==2;0])>0) find(diff([0;estimateLabel(:,1)==2;0])<0)];
        kp = diff(sz_ep,[],2)>10;
        sz_ep = sz_ep(kp,:);
        
        save(fout,'estimateLabel','ops','ClassifierFileOutput','sz_ep')
   % end
end


%%
plot(estimateLabel)
hold on
plot(estimateLabel(:,1))