%% Load Data
load('R:\DGregg\NeuralData\seizurePaths_annotatedMovement.mat');
%% Set Parameters

bp        = [.1 300];
preTime   = 120;   % s of pre-ictal LFP loaded
postTime  = 600;   % s of post-ictal LFP loaded
fixedWinS = 120;   % s: baseline estimation window at the start of preTime
dt  = 10;    % s: bin width for post-ictal power traces
nPostBins = postTime/dt;
outlierMAD = 5;    % samples beyond this many MAD-equivalent SDs from the
                    % preictal median are treated as spikes and excluded
                    % from the baseline mean

nCh = 8;
chName = {'M2 (L)','M2 (R)','Thal (L)','BLA (R)', ...
          'CA1 (L)','CA1 (R)','LDT (L1)','LDT (L2)'};
%% data extraction

clear pctSuppBin spont
for i = 1:length(seizure)
    % i
    lfpfil = seizure(i).info.sessiondata.lfp_file;
    xmlfil = strrep(lfpfil,'lfp','xml');
    xml    = LoadXml(xmlfil);
    fs     = xml.lfpSampleRate;
    nChannelsXML  = xml.nChannels;
    ChannelsSubj  = seizure(i).info.sessiondata.channelID;

    fixedWin = 1:(fixedWinS*fs);
    

    for j = 1:length(seizure(i).data)
        spont{i}(j) = strcmp(char(seizure(i).data(j).spontaneous),'spontaneous');

        szstart = seizure(i).data(j).szstart;
        szdur_j = seizure(i).data(j).szdur;
        idx_sz_end = ceil(fs*(preTime+szdur_j));

        for k = 1:length(ChannelsSubj)
            lfp = LoadBinary(lfpfil,'frequency',fs, ...
                'nChannels',nChannelsXML,'start', szstart-preTime, ...
                'duration',preTime+szdur_j+postTime, ...
                'channels',ChannelsSubj(k));
            lfp = double(lfp);
            IA  = InstAmplitude(BandpassFilter(lfp,fs,bp));

            IA_pre = IA(1:preTime*fs);

            w      = IA_pre(fixedWin);
            medW   = nanmedian(w);
            madW   = mad(w(~isnan(w)),1);              % median absolute deviation
            isSpike = abs(w - medW) > outlierMAD*1.4826*madW;
            bl = nanmean(w(~isSpike));
            if isnan(bl) || all(isSpike)
                bl = medW;  
                disp('flagged everything as outlier for baseline\n')
            end

            % NOTE: This is a change to suppression instead of PI power!!!!! 
            pctSupp = -(IA - bl)/bl;
            pctSupp = pctSupp(idx_sz_end:end);
           % pctSupp = pctSupp(~isnan(pctSupp));
            pctSupp = nanPad(pctSupp', postTime*fs + 1);  

            ts = (1:length(pctSupp))/fs;
            pctSuppBin{i}(j,:,k) = avghist(ts, pctSupp, 1:dt:ts(end));
        end
    end
    % i
    fprintf('extracted session %d/%d\n', i, length(seizure));
end


%% subject grouping to find IHKA
clear subj
for i = 1:length(seizure)
    subj{i} = seizure(i).subject{1};
end
IHKA = contains(subj,'EDS');


%% Pull out censoring data
for i = 1:length(seizure)
    evokedR = false;   % tracks "in an evoked run" while scanning forward
    evokedL = false;   % same, while scanning backward is handled inline below
    if isempty(seizure(i).data), continue, end

    for j = 1:length(seizure(i).data)
        seizure(i).data(j).pctSupp = pctSuppBin{i}(j,:,:);
        seizure(i).data(j).pctSupp_mean = nanmean(pctSuppBin{i}(j,1:end-20,6),2);
        % seizure(i).data(j).excitability = excitabilityBin{i}(j,:);

        spontStr = char(seizure(i).data(j).spontaneous);
        
        if strcmp(spontStr,'spontaneous') && ~evokedR && j < length(seizure(i).data) ...
                && strcmp(char(seizure(i).data(j+1).spontaneous),'spontaneous') ...
                && strcmp(seizure(i).data(j).censor.right_type,'uncensored')
            seizure(i).data(j).next_seizure_interval = seizure(i).data(j+1).szstart - seizure(i).data(j).szstart;
            seizure(i).data(j).censoredR = false;
        elseif strcmp(spontStr,'spontaneous') && ~evokedR && j < length(seizure(i).data) ...
                && ~strcmp(seizure(i).data(j).censor.right_type,'uncensored')
            seizure(i).data(j).next_seizure_interval = seizure(i).data(j).censor.right_time - seizure(i).data(j).szstart;
            seizure(i).data(j).censoredR = true;
        elseif strcmp(spontStr,'spontaneous') && ~evokedR && j == length(seizure(i).data)
            seizure(i).data(j).next_seizure_interval = seizure(i).data(j).censor.right_time - seizure(i).data(j).szstart;
            seizure(i).data(j).censoredR = true;
        elseif strcmp(spontStr,'evoked')
            seizure(i).data(j).next_seizure_interval = nan;
            seizure(i).data(j).censoredR = true;
            evokedR = true;
        else
            seizure(i).data(j).next_seizure_interval = nan;
            seizure(i).data(j).censoredR = true;
        end

       
        if strcmp(spontStr,'spontaneous') && ~evokedL && j > 1 ...
                && strcmp(char(seizure(i).data(j-1).spontaneous),'spontaneous') ...
                && strcmp(seizure(i).data(j).censor.left_type,'uncensored')
            seizure(i).data(j).last_seizure_interval = seizure(i).data(j).szstart - seizure(i).data(j-1).szstart;
            seizure(i).data(j).censoredL = false;
        elseif strcmp(spontStr,'spontaneous') && ~evokedL && j > 1 ...
                && ~strcmp(seizure(i).data(j).censor.left_type,'uncensored')
            seizure(i).data(j).last_seizure_interval = seizure(i).data(j).szstart - seizure(i).data(j).censor.left_time;
            seizure(i).data(j).censoredL = true;
        elseif strcmp(spontStr,'spontaneous') && ~evokedL && j == 1
            seizure(i).data(j).last_seizure_interval = seizure(i).data(j).szstart - seizure(i).data(j).censor.left_time;
            seizure(i).data(j).censoredL = true;
        elseif strcmp(spontStr,'evoked')
            seizure(i).data(j).last_seizure_interval = nan;
            seizure(i).data(j).censoredL = true;
            evokedL = true;
        else
            seizure(i).data(j).last_seizure_interval = nan;
            seizure(i).data(j).censoredL = true;
        end
    end
end


%% build analysis table and fit the model
% onset:   1 = LVF, 2 = HYP
% behavioral: 1 = convulsive, 2 = not convulsive, 3 = failed
% movement:   1 = moving <4s pre-sz, 2 = not moving, 3 = failed

seizure1 = seizure(IHKA);
onsetType = []; 
behavioral = [];
movement = []; 
ses_ID = [];
is_evoked = [];
isSpont = []; 
generalized = [];
P = [];   

% Model: y = A * exp(-t/tau)
ft = fittype('A*exp(-x/tau)+c', ...
    'independent','x', ...
    'coefficients',{'A','tau','c'});


s_id = [];% (session, seizure)

for i = 1:length(seizure1)
    if isempty(seizure1(i).data), continue, end
    for j = 1:length(seizure1(i).data)
            switch seizure1(i).data(j).behavioral
                case categorical({'behavioral'})
                    behavioral(end+1,1) = 1;
                case categorical({'non-behavioral'})
                    behavioral(end+1,1) = 2;
                otherwise
                    behavioral(end+1,1) = 3;
            end
            switch seizure1(i).data(j).movement
                case categorical({'movement'})
                    movement(end+1,1) = 1;
                case categorical({'non-movement'})
                    movement(end+1,1) = 2;
                case categorical({'failed'})
                    movement(end+1,1) = 3;
            end
            switch seizure1(i).data(j).onsetType
                case categorical({'LVF'})
                    onsetType(end+1,1) = 1;
                case categorical({'HYP'})
                    onsetType(end+1,1) = 2;
            end
            switch seizure1(i).data(j).spontaneous
                case categorical({'spontaneous'})
                    is_evoked(end+1,1)  = 1;
                case categorical({'evoked'})
                    is_evoked(end+1,1)  = 2;
                otherwise
                    is_evoked(end+1,1)  = 3;
            end
            ses_ID = [ses_ID; i j]; 
            generalized = [generalized; seizure1(i).data(j).generalized == 1];
            isSpont(end+1,1) = strcmp(char(seizure1(i).data(j).spontaneous),'spontaneous');
            P(end+1,:,:) = seizure1(i).data(j).pctSupp;
            t = [(1:size(P,2))*10]';
            s_id(end+1) = i;


            for k = 1:size(P,3)
                PI_power = -seizure1(i).data(j).pctSupp(1,:,k)';
                y = PI_power;


                % Initial guesses
                A0 = y(1);
                tau0 = (max(t)-min(t))/3;
                c0 = 0;

                opts = fitoptions(ft);
                opts.StartPoint = [A0 tau0 c0];

                % A can be positive or negative; tau must be positive.
                % NOTE: y is on the percent-suppression scale (pctSupp is
                % -100*(IA-bl)/bl), not the old raw-fraction scale, so the
                % offset bound is widened accordingly (was +/-10, sized for
                % a +/-1-ish fractional y).
                opts.Lower = [-Inf 0 -100];
                opts.Upper = [ Inf 100 100];
                kp_nan = ~isnan(y);
                [f,gof] = fit(t(kp_nan),y(kp_nan),ft,opts);


                seizure1(i).data(j).A(k)   = f.A;
                seizure1(i).data(j).tau(k) = f.tau;
                seizure1(i).data(j).c(k) = f.c;
                seizure1(i).data(j).R2(k) = gof.rsquare;
                seizure1(i).data(j).Y_hat(:,k) = f.A.*exp(-t/f.tau)+f.c;
            end
    end
end
%%

% Pulled from seizure1, not seizure(IHKA): seizure1 = seizure(IHKA)
censoredR = cell2mat(arrayfun(@(a) cell2mat({a.data.censoredR}'), seizure1, 'uni',0)');
censoredL = cell2mat(arrayfun(@(a) cell2mat({a.data.censoredL}'), seizure1, 'uni',0)');
timeNext  = cell2mat(arrayfun(@(a) cell2mat({a.data.next_seizure_interval}'), seizure1, 'uni',0)');
timeLast  = cell2mat(arrayfun(@(a) cell2mat({a.data.last_seizure_interval}'), seizure1, 'uni',0)');
szdur     = cell2mat(arrayfun(@(a) cell2mat({a.data.szdur}'), seizure1, 'uni',0)');
A_fit     = cell2mat(arrayfun(@(a) cell2mat({a.data.A}'), seizure1, 'uni',0)');
tau_fit   = cell2mat(arrayfun(@(a) cell2mat({a.data.tau}'), seizure1, 'uni',0)');
c_fit     = cell2mat(arrayfun(@(a) cell2mat({a.data.c}'), seizure1, 'uni',0)');
R2_fit    = cell2mat(arrayfun(@(a) cell2mat({a.data.R2}'), seizure1, 'uni',0)');

onsetHYP = onsetType == 2;

%% 

clear Beta p_val
for b = 1:2
for k = 1:8
kpN = ~isnan(timeNext) & behavioral==b &is_evoked==1 ; % analysze time to next
X = [c_fit(kpN,k)  szdur(kpN) ];
Y = timeNext(kpN);

censored = [censoredR(kpN)];

coxModel = fitcox(X, Y, 'Censoring', censored);
Beta(k,b,:)  =coxModel.Coefficients.Beta;
p_val(k,b,:) =       coxModel.Coefficients.pValue;
end
end
close all

       figure 
bar(Beta(:,:,1))
hold on
sig1 = p_val(:,1,1);
sig1(sig1>.05) = nan;
sig1(~isnan(sig1)) = 1;

sig1a = nan(8,1);
sig1a(p_val(:,1,1)<.1 & p_val(:,1,1)> .05) = 1;




sig2 = p_val(:,2,1);
sig2(sig2>.1 ) = nan;
sig2(~isnan(sig2)) = 1;

text(1.15:8.15,-sig2*2,'#','color','r')
plot(.85:7.85,sig1*2,'*','color','b')
text(.85:7.85,sig1a*2,'#','color','b')
legend({'Conv.','Non-Conv.','',''})
ylim([-2.5 2.5])
ylabel('offset')
set(gca,'xtick',1:8,'xticklabel',chName)


%%
clear beta_cox p_cox


for t = 1:nPostBins
    for ch = 1:8
        tmp = [];

        for  i = 1:(length(seizure1))

            if ~isempty(seizure1(i).data)
                for  j = 1:length(seizure1(i).data)

                    tmp = [tmp;seizure1(i).data(j).pctSupp(:,t,ch)];

                end
            end
        end
        x = [tmp] ;

        if ch ==1
            x1 = x;
        end

        kpN = ~isnan(timeNext); % analysze time to next
        kpL = ~isnan(timeLast);% analysze time from past
      
        % now also drop rows where x is NaN, not just where time is NaN
        kpN = kpN & ~isnan(x);
        kpL = kpL & ~isnan(x);
       

        %tbl = table(timeNext(kp2), censored(kp2), x(kp2),timeLast(kp2),szdur(kp2), 'VariableNames', {'time','censored','PI','timeLast','szdur'});
        %coxModel = fitcox(x(kp2), timeNext(kp2), 'Censoring', censored(kp2));
        
%         X = [x(kpL) szdur(kpL)];
%         Y = timeLast(kpL);
%         censored = [censoredL(kpL)];
       
         
         
          kp_all = kpN ;
          behavioral1= behavioral==1;
        %  X = [zscore(x(kp_all)) szdur(kp_all) ];
         X = [zscore(x(kp_all)) behavioral1(kp_all) behavioral1(kp_all).*zscore(x(kp_all)) szdur(kp_all) ];
        Y = timeNext(kp_all);
       
        censored = [censoredR(kp_all)];

        coxModel = fitcox(X, Y, 'Censoring', censored);
        beta_cox(ch,t,:) = coxModel.Coefficients.Beta;
        p_cox(ch,t,:) = coxModel.Coefficients.pValue;
    end

end


%%
close all
figure
ts_plot = 1:10:ts(end);
ax  = tight_subplot(4,2);
for i = 1:8
    
    ix = find(p_cox(i,:,3)<.05);
    axes(ax(i))
    plot(ts_plot,exp(beta_cox(i,:,3)))
    hold on
    if any(ix)
        plot(ts_plot(ix),1.8,'.','markersize',10,'color','k')
    end
    xlim([0 300])
    ylim([.25 5])
    
    title(chName{i})
    if ~(i==7 | i==8)
        set(gca,'xticklabel','')
    else
        xlabel('Time from seizure end (s)')
    end
    
    ylabel('Hazard ratio')
end



%% PI timecourse
colLVF = [100/255,170/255,215/255];
colHYP = [233/255,152/255,117/255];

ts_PETH_bin = (1:nPostBins)*dt;
kpLVF = is_evoked==1 & (behavioral==1);
kpHYP = is_evoked==1 & behavioral==2;

% fprintf('\npostictal timecourse: n spontaneous Convulsive = %d, non-Convulsive = %d\n', ...
%     sum(kpLVF), sum(kpHYP));

figure
ax = tight_subplot(4,2);
for ch = 1:nCh
    axes(ax(ch)) 
    plotMeanSEM(ts_PETH_bin, -P(kpLVF,:,ch)*100, colLVF);

    hold on
    plotMeanSEM(ts_PETH_bin, -P(kpHYP,:,ch)*100, colHYP);

    
    clear pp
    for h = 1:nPostBins
        a = P(kpLVF,h,ch); a = a(~isnan(a));
        b = P(kpHYP,h,ch); b = b(~isnan(b));
        if numel(a) > 1 && numel(b) > 1
            [~,pp(h)] = ttest2(a,b);
        else
            pp(h) = nan;
        end
    end
    pp(pp>.01) = nan;
    pp(~isnan(pp)) = 60;
    plot(ts_PETH_bin,pp,'k','linewidth',6)
    xlim([0 120])
    title(chName{ch})
    % xlim([0 300])
    yline(0)
    ylim([-70 70])
    set(gca,'fontsize',16)
    if ch==1
        legend({'Convulsive','','Non-Convulsive'},'Location','best')
    end
    if ~(ch==7 || ch==8)
        set(gca,'xticklabel','')
    end
end

%% Fit Cox model by behavior
intWin = 1:nPostBins;   
conditions = struct( ...
    'name', {'Convulsive','non-Convulsive'}, ...
    'mask', { ...
    behavioral==1, ...
    behavioral==2 ...
    % generalized==1, ...
    % movement==1, ...
    % movement==2 ...
    });

clear beta_ctx se_ctx p_ctx n_ctx ev_ctx
for ch = 1:nCh
    x = nan(size(timeNext));
    row = 0;
    for i = 1:length(seizure)
        if isempty(seizure(i).data) || ~IHKA(i), continue, end
        for j = 1:length(seizure(i).data)
            row = row + 1;
            x(row) = nanmean(seizure(i).data(j).pctSupp(:,intWin,ch),2);
        end
    end

    for c = 1:numel(conditions)
        kp = conditions(c).mask(:) & ~isnan(timeNext(:)) & ~isnan(x(:)) & ~isnan(szdur(:));

        n_ctx(c,ch)  = sum(kp); 
        ev_ctx(c,ch) = sum(~censoredR(kp)); 
        if n_ctx(c,ch) < 10 || ev_ctx(c,ch) < 5
            beta_ctx(c,ch) = nan; se_ctx(c,ch) = nan; p_ctx(c,ch) = nan;
            continue
        end

        xz = zscore(x(kp));
        dz = zscore(szdur(kp));
        X  = [xz, dz];

        m = fitcox(X, timeNext(kp), 'Censoring', censoredR(kp));
        beta_ctx(c,ch) = m.Coefficients.Beta(1);   % suppression term
        se_ctx(c,ch)   = m.Coefficients.SE(1);    
        p_ctx(c,ch)    = m.Coefficients.pValue(1); 
    end
end

%% Plot cox by onset type 

figure
HR   = exp(beta_ctx);                            % nConditions x nCh
CIlo = exp(beta_ctx - 1.96*se_ctx);
CIhi = exp(beta_ctx + 1.96*se_ctx);

condColors = [colLVF; colHYP];   % extend if you add more conditions above

b = bar(1:nCh, HR', 'grouped');
for c = 1:numel(conditions)
    b(c).FaceColor = condColors(c,:);
end
hold on

nbars = numel(conditions);
groupwidth = min(0.8, nbars/(nbars+1.5));
for c = 1:numel(conditions)
    xpos = (1:nCh) - groupwidth/2 + (2*c-1)*groupwidth/(2*nbars);
    errorbar(xpos, HR(c,:), HR(c,:)-CIlo(c,:), CIhi(c,:)-HR(c,:), ...
        'k','linestyle','none','linewidth',1.5)
    sig = find(p_ctx(c,:) < .05);
    if any(sig)
        plot(xpos(sig), CIhi(c,sig)+0.1, '*k','markersize',10)
    end
end
yline(1,'--k')

set(gca,'xtick',1:nCh,'xticklabel',chName,'fontsize',12)
xtickangle(45)
ylabel('Hazard ratio (per SD postictal suppression)')
legend({conditions.name},'Location','best')
title('Suppression''s effect on next-seizure risk, fit separately by behavioral manifestation')


