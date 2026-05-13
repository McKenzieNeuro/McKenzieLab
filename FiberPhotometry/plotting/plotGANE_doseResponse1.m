%%

dirs = [...
    {'R:\DANEHippocampalResponse\NE2h23\DoseResponse\NE2h23_250903'};...
    {'R:\DANEHippocampalResponse\NE2h22\DoseResponse\NE2h22_250903'};...
    {'R:\DANEHippocampalResponse\NE2h23\DoseResponse\NE2h23_250902'};...
    {'R:\DANEHippocampalResponse\NE2h22\DoseResponse\NE2h22_250902'};...
    {'R:\DANEHippocampalResponse\NE2h24\doseResponse\NE2h24_250902'};...
    {'R:\DANEHippocampalResponse\NE2h24\DoseResponse\NE2h24_250903'};...
    %{'R:\DANEHippocampalResponse\NE2h22\AM251\NE2h22_250903'};...
    ];



%%
%sig_mean =[];pulseID=[];
for i = 1:length(dirs)
    
    cd(dirs{i})
    [~,subj] = fileparts(fileparts(fileparts(dirs{i})));

    switch subj
        case 'NE2h22'
            ch = 8;
        case 'NE2h23'
            ch = 7;
        case 'NE2h24'
             ch = 7;
    end

    %find TDT directory
    dirN = fileparts(getAllExtFiles(dirs{i},'Tbk',1));
    %data = TDTbin2mat(dirN);
    [signal_DFoF2,ts_data2,fs,data] = sm_getSignal_DFoF(dirN{1});
    k = gaussian2Dfilter([fs 1],fs);
    signal_DFoF = nanconvn(signal_DFoF2,k');
    % signal_DFoF = nanconvn(data.streams.x405A.data,k');
    ts_data2 = (1:length(data.streams.x405A.data))/fs;
    %get pulses
    digF = getAllExtFiles(dirs{i},'dat',1);
    digF = digF(contains(digF,'digitalin'));
    pulse_fil = [fileparts(digF{1}) filesep 'TTL_pulse.mat'];
    
    if ~exist(pulse_fil)
        [ups,dwns]  = sm_getDigitalin(fileparts(digF),'digitalin.dat',30000,16);
    else
        load(pulse_fil)
    end
    TDT_on = data.epocs.Pu1_.onset;
    
    %convert TDT time into Intan time
    
    Intan_ts_FP = interp1(TDT_on,dwns{3}(1:length(TDT_on)),ts_data2,'linear','extrap');
    
    %get pulse times
    fname = strrep(digF{1},'digitalin','amplifier_analogin_auxiliary_int16');
    pulse_fil = [fileparts(digF{1}) filesep 'pulseInfo.mat'];
    if ~exist(pulse_fil)
        [onset,offset,pulse,center] = sm_getPulseTime(fname,9,'pulse');
        
        % do some QC
        
        % should be >5s apart and ~10ms
        
        kp1 = [true;diff(onset{1}(:,1))>.05];
        % kp1 = [0;diff(onset{1,1}(:,1))]>4; % pulses should be more than 4s apart
        kp2 = (offset{1} - onset{1}) >.08 & (offset{1} - onset{1}) <.62;
        
        kp = find(kp1& kp2);
        if any (kp)
            %save and categorize by pulse amplitude
            pulseInfo(1).time = [onset{1}(kp) offset{1}(kp)];
            pulseInfo(1).amp = pulse{1}(kp);
            
            % number of intensities, Change it before running!!!!!!!!!!
            
            
            ix = kmeans(pulseInfo(1).amp',numpulses);
            pulse_amp = accumarray(ix,pulseInfo(1).amp',[numpulses 1],@nanmean,nan);
            [~,b]  =sort(pulse_amp);
            %reassign
            pulseInfo(1).id = ix;
            for ii = 1:numpulses
                pulseInfo(1).id(ix==b(ii)) = ii;
            end
            
            save([fileparts(fname) filesep 'pulseInfo.mat'],'pulseInfo')
            
            
        end
    else
        load(pulse_fil)
    end
    
    pulseTime_FP = interp1(Intan_ts_FP,ts_data2,pulseInfo.time(:,1),'linear','extrap');
    
    [ix,early,late,ts] = sm_getIndicesAroundEvent(pulseTime_FP,20,120,fs,length(signal_DFoF2));
    kp = ~early & ~late;
  

    [~,b0] = bestmatch(-10,ts);
    [~,b1] = bestmatch(0,ts);
    [~,b2] = bestmatch(20,ts);
    
    sig_mean{i} = signal_DFoF(ix(kp,:));

    AUC_FP{i} = nanmean(sig_mean{i}(:,b1:b2),2) -  nanmean(sig_mean{i}(:,b0:b1),2);


    pulseID{i} =  pulseInfo(1).amp(kp);

    lfp  = [fileparts(digF{1}) filesep 'amplifier_analogin_auxiliary_int16.dat'];

    ts_STIM  =pulseInfo.time(kp,1);
    lfp_pulse{i} = [];
    for j = 1:length(ts_STIM)
        lfp_pulse{i}(j,:) = LoadBinary(lfp,'nchannels',12,'channels',ch,'frequency',30000,'start',ts_STIM(j)-.1,'duration',.6);
        baseline = mean( lfp_pulse{i}(j,1:2900));
        peak = mean(lfp_pulse{i}(j,3301:6000));
        
        

        lfp_slope{i}(j) = (.195*(peak - baseline))/ .02;

    end
end

ts_lfp  =(1:length(  lfp_pulse{1}))/30000-.1;


%%

id = cell2mat(pulseID)';

ix = kmeans(id,9);
pulse_amp = accumarray(ix,id,[9 1],@nanmean,nan);
[a,b]  =sort(pulse_amp);
%reassign
id1 = ix;
for ii = 1:9
    id1(ix==b(ii)) = ii;
end
%%
numpulses = 9;
           
           all_resp = cell2mat(sig_mean');
            all_resp_lfp = cell2mat(lfp_pulse');
           gr = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),sig_mean',num2cell([1:length(sig_mean)]',2),'uni',0));
           
%close all
col = linspecer(numpulses,'jet');
figure
for i = 1:9
    
    plot(ts,nanmean(all_resp(id1==i & gr~=1,:)),'color',col{i})
    hold on
end
colormap(cell2mat(col))
colorbar

%%
[~,b0] = bestmatch(-10,ts);
[~,b1] = bestmatch(0,ts);
[~,b2] = bestmatch(5,ts);
[~,b3] = bestmatch(60,ts);
[~,b4] = bestmatch(90,ts);


[~,b0a] = bestmatch(-.1,ts_lfp);
[~,b1a] = bestmatch(0,ts_lfp);
[~,b2a] = bestmatch(.002,ts_lfp);
[~,b3a] = bestmatch(.025,ts_lfp);



clear u_rep_early u_rep_vearly u_rep_late u_rep_vlate
for i = 1:numpulses
    
   
    u_rep_vearly{i} = nanmean(all_resp(id1==i & gr~=7,b0:b1),2);
    u_rep_early{i} = nanmean(all_resp(id1==i& gr~=7,b1:b2),2);
    u_rep_late{i} = nanmean(all_resp(id1==i& gr~=7,b2:b3),2);
    u_rep_vlate{i} = nanmean(all_resp(id1==i& gr~=7,b3:b4),2);
    u_rep_vearly_am251{i} = nanmean(all_resp(id1==i & gr==7,b0:b1),2);
    u_rep_early_am251{i} = nanmean(all_resp(id1==i& gr==7,b1:b2),2);


    u_rep_vearly_lfp{i} = nanmean(all_resp_lfp(id1==i & gr~=7,b0a:b1a),2);
    u_rep_early_lfp{i} = nanmean(all_resp_lfp(id1==i& gr~=7,b2a:b3a),2);



    u_FP_AUC{i} =  nanmean(all_resp(id1==i& gr~=7,b1:b2),2) -  nanmean(all_resp(id1==i & gr~=7,b0:b1),2);
    u_LFP_slope{i} =  nanmean(all_resp_lfp(id1==i& gr~=7,b2a:b3a),2) -  nanmean(all_resp_lfp(id1==i & gr~=7,b0a:b1a),2);

end

figure
X_u = cellfun(@median,u_LFP_slope(3:end));
Y_u = cellfun(@median,u_FP_AUC(3:end));

X_SEM= cellfun(@SEM,u_LFP_slope(3:end));
Y_SEM = cellfun(@SEM,u_FP_AUC(3:end));

H=errorbarxy(X_u(1:end),Y_u(1:end),X_SEM(1:end),Y_SEM(1:end),{'ko-', 'k', 'k'});

xlabel('mean LFP - 0-25 ms')
ylabel('mean FP - 0-5 s')


hold on
plot(cell2mat(u_LFP_slope'),cell2mat(u_FP_AUC'),'o')
set(gca,'fontsize',16)
%%
close all

ds = 100;
figure

 
shg
hold on

plotMeanSEM(ts_lfp(1:ds:end),.195*(all_resp_lfp(id1==3& gr~=7,(1:ds:end))),'b')
plotMeanSEM(ts_lfp(1:ds:end),.195*(all_resp_lfp(id1==5& gr~=7,(1:ds:end))),'g')
plotMeanSEM(ts_lfp(1:ds:end),.195*(all_resp_lfp(id1==9& gr~=7,(1:ds:end))),'r')
set(gca,'fontsize',16)
figure

 
shg
hold on

plotMeanSEM(ts(1:ds:50000),(all_resp(id1==3& gr~=7,(1:ds:50000))),'b')
plotMeanSEM(ts(1:ds:50000),(all_resp(id1==5& gr~=7,(1:ds:50000))),'g')
plotMeanSEM(ts(1:ds:50000),(all_resp(id1==9& gr~=7,(1:ds:50000))),'r')
set(gca,'fontsize',16)

%%
col = linspecer(5,'jet');
figure
for i = 5:-1:1
    
    plot(ts,nanmean(sig_mean(i,:,:),3),'color',col{i},'linewidth',3)
    hold on
end

colormap(cell2mat(col))
colorbar
%%

col = linspecer(5,'jet');
figure
for i = 1:5
    
    plot(ts,sig_mean465(i,:),'color',col{i})
    hold on
end


