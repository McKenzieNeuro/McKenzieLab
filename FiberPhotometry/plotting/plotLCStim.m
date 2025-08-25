%%

%fixed stim

dirs = [ ...
    {'R:\DANEHippocampalResponse\NE2h18\Opto Stim\NE2h18_250623'};...
    {'R:\DANEHippocampalResponse\NE2h19\Opto Stim\NE2h19_250624'};...
    {'R:\DANEHippocampalResponse\NE2h18\Opto Stim\NE2h18_250627'};...
    {'R:\DANEHippocampalResponse\NE2h20\Opto Stim\NE2h20_250627'};...
    {'R:\DANEHippocampalResponse\NE2h20\Opto Stim\NE2h20_250630'};...
    
    ];

%%
for i = 1:length(dirs)
    
    cd(dirs{i})
    
    %get mouse
    idx = regexp(dirs{i},'NE2h');
    mouse = dirs{i}(idx(1):idx(1)+5);
    
    switch mouse
        case 'NE2h18'
            ch = 29;
        case 'NE2h19'
            ch = 22;
        case 'NE2h20'
            ch = 23;
    end
    
    % find TDT directory
    
    
   
     fils=  getAllExtFiles(dirs{i},'Tbk',1);
    TDTdir = fileparts(fils{1});
   
    [signal_DFoF2,ts_data2,fs] = sm_getSignal_DFoF(TDTdir);
    data = TDTbin2mat(TDTdir);
    %load stims
    stim  = data.epocs.Pe1_.onset;
    
    
    % find intan dir
    fils=  getAllExtFiles(pwd,'dat',1);
    Intandir = fileparts(fils{1});
    
    TTLfil = [Intandir filesep 'TTL_pulse.mat'];
    
    if exist(TTLfil)
        load(TTLfil)
        
    else
        
        [ups,dwns]  = sm_getDigitalin(Intandir,'digitalin.dat',30000,16);
    end
    
    stim_intan = ups{1}(diff([0;ups{1}])>1);
    
    % get evoked
    
    
    
    [ix,early,late,ts] = sm_getIndicesAroundEvent(stim,100,200,fs,length(signal_DFoF2));
    [~,b1] = bestmatch(0,ts);
    [~,b2] = bestmatch(10,ts);
    kp = ~early & ~late;
    stim_intan = stim_intan(kp);
    ix = ix(kp,:);
    
    
    stim = stim(kp);
    
    k = gaussian2Dfilter([fs 1],fs);
    signal_DFoF = nanconvn(signal_DFoF2,k');
    
    NE_peth = signal_DFoF(ix);
    
    % fit
    
    x=ts(b2:end);
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
        'Lower', [-inf], 'Upper', [ 0]);
    
    fitfun = 'exp(a*x)';
    
    
    
    clear params NE_peth1 pk base
    for ii = 1:length(stim)
        y  = NE_peth(ii,:);
        
        base(ii) = nanmean(y(1:b1));
        y = y(b2:end)-nanmean(y(1:b1));
        pk(ii)= nanmean(y(1:10));
        y = y/nanmean(y(1:10));
        
        NE_peth1(ii,:) = y;
        [f,gof] = fit(x(:),y(:),fitfun,options);
        params(ii,:) = coeffvalues(f);
        rmse(ii) = gof.rmse;
    end
    %get spectra
    
    kp =params>median(params);
    lfpfil = [Intandir filesep 'amplifier_analogin_auxiliary_int16.lfp'];
  
    if ~exist(lfpfil)
          bz_LFPfromDat(Intandir,'basename','amplifier_analogin_auxiliary_int16')
    end
    figure
    [h,wavspec,ts,freqs] = sm_EventTriggeredSpect(lfpfil,stim_intan,'channel',ch,'freqs',logspace(log10(2),log10(300),50));
    
    
    figure
    [~,b] = histc(params,-.02:.002:-.006);
    col = linspecer(max(b),'jet');
    
    base_lfp = nanmedian(wavspec(:,1000:62500,:),2);
    base_lfp = squeeze(base_lfp);
    for ii  = 1:max(b)
        
        semilogx(freqs,nanmean(base_lfp(:,b==ii&rmse'<.15),2),'color',col{ii},'linewidth',2)
        hold on
    end
    
    colormap(cell2mat(col))
    cbh = colorbar; % Get the handle to the current colorbar
    
    numTicks = 8; % Set the desired number of ticks
    cbh.Ticks = linspace(0, 1, numTicks); % Create tick locations from 0 to 1
    cbh.TickLabels = num2cell(-.05:.005:0);
    xlabel('frequency')
    ylabel('power')
    cbh.Label.String = 'NE dedcay';
    
    
    %        plot(f,x,y)
    %     gof.rmse
    %        waitforbuttonpress
    %        close all
    
    close all
    col = flipud(linspecer(length(stim),'jet'));
    [~,b] = sort(params(:,1));
    figure
    for ii = length(stim):-1:1
        if rmse(ii)<.25
            plot(x,NE_peth1(b(ii),:),'color',col{ii})
            
        end
        hold on
    end
    
    col = linspecer(length(stim),'jet');
    [~,b] = sort(params(:,1));
    figure
    for ii = 1:length(stim)
        if rmse(ii)<.25
            plot(x,NE_peth1(b(ii),:),'color',col{ii})
            
        end
        hold on
    end
end




%%
% load data


%%


%%



%%


%%

base_lfp = nanmedian(wavspec(:,1000:62500,:),2);
base_lfp = squeeze(base_lfp);

d = tsne(base_lfp');
d1= sortby(d,params);
close all
figure
col = linspecer(33,'jet');
for i  = 1:33
    plot(d1(i,1),d1(i,2),'.','markersize',100,'color',col{i})
    hold on
    
end

%%

ix = kmeans(d,2);
figure
plot(ts,nanmean(NE_peth(ix==1,:)))
hold on
plot(ts,nanmean(NE_peth(ix==2,:)))

%%
figure
plot(ts,nanmean(NE_peth(params(:,1)<-.015 &rmse'<.25,:)))
hold on
plot(ts,nanmean(NE_peth(params(:,1)>-.015 &rmse'<.25,:)))

%%
figure
kp =params>median(params);
semilogx(freqs,nanmean(base_lfp(:,~kp),2))
hold on
plot(freqs,nanmean(base_lfp(:,kp),2),'k')

%%


