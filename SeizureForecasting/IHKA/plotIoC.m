% "initialize" variables that will load data across sessions in a loop
%F1-Scoring, sensitivity, specificity,
%area under the receiver‐operating characteristic curve and time in false warning
ID_IoC =[];C =[];IoC =[];p_sz =[];BS =[];RSS = [];
TiW = []; TiW_r = [];seizure_ITI = []; AUC =[]; f1 =[];
sensitivity = [];specificity =[];
idx= 1 ;
all_pred = [];
all_tim2seiz = [];
all_IDs = [];
prev = [];
precision = [];
predictedR = [];
predicted = [];

alarmsPerHour = [];
sens_block = [];
TiW_alarmOnly = [];
blockLenMean = [];
chance_sens = [];
AUC1 =[];
AUC2= [];
N1 = [];
N2 = [];
nBlocks =[];
%loop over all predictions where each 1-s time window has been categorized
%by which time window it is from

nGr = 6;
clear ok1 d
for j = 1:72
    
    %load predictions for each 24-hr period
    load(['R:\IHKA_gross\prediction\predict_' num2str(j) '.mat']);
    
    % the true label can be shorter that the estimate since the esimtate
    % requires 2-s of data, so yolk them to the same length
    trueLabel = trueLabel(1:length(estimateLabel));
    
    %build time stamps
    ts = 1:length(time2seizure);
    Est = nan(nGr,length(estimateLabel));
    for i = 1:nGr
        
        Est(i,:) = estimateLabel(:,1)==i;
        
    end
    [~,ix] =  histc(seizure_start,ts);
    
    ix1 = repmat(ix,1,7201)+repmat(-3600:3600,length(ix),1);
    kp = all(ix1>0 & ix1<=size(Est,2),2);
    ix1 = ix1(kp,:);
    ix = ix(kp);
    tims = seizure_start(kp);
    for i = 1:length(ix)
        if i ==1
            actual_Y = trueLabel(1:ix(i)+600);
            predicted_Y = estimateLabel(1:ix(i)+600);
            confLabel_Y = estimateLabel(1:ix(i)+600,2:end);
        else
            actual_Y = trueLabel(ix(i-1)+600:ix(i)+600);
            predicted_Y = estimateLabel(ix(i-1)+600:ix(i)+600);
            confLabel_Y = estimateLabel(ix(i-1)+600:ix(i)+600,2:end);
        end
        predicted_Y = predicted_Y(actual_Y>0);
        confLabel_Y = confLabel_Y(actual_Y>0,:);
        actual_Y = actual_Y(actual_Y>0);
        
        ixx  =unique([actual_Y(:);predicted_Y(:)]);
        tmp1 = nan(nGr,nGr);
        
        tmp1(ixx,ixx) = confusionmat(actual_Y,predicted_Y);
        
        
        for jj = 1:6
            if any(actual_Y==jj) && ~all(actual_Y==jj)
                [X,Y,T,AUCt(jj)] = perfcurve(actual_Y==jj,confLabel_Y(:,jj),1);
            else
                AUCt(jj) = nan;
            end
            [f1t(jj), sensitivityt(jj), specificityt(jj),precisiont(jj)] = computeClassificationMetrics(actual_Y==jj,predicted_Y==jj);
            prevt(jj) = mean(actual_Y==jj);   % inside the seizure loop
        end
        f1 = [f1;f1t];
        sensitivity = [sensitivity;sensitivityt];
        specificity = [specificity;specificityt];
        precision = [precision;precisiont];
        AUC = [AUC;AUCt];
        prev = [prev;prevt];
        % we have an observed confusion matrix
        C = cat(3,C,tmp1);
        nBoot = 1000;
        Cr = nan(nGr,nGr,nBoot);BSr = nan(nBoot,nGr);
        for ii = 1:nBoot
            
            % get conf
            actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
            tmp = confusionmat(actual_Yr,predicted_Y);
            
            tmp = tmp./sum(tmp,2);
            
            predictedRt = any(actual_Y(actual_Yr==4)==4);
            
            Cr(ixx,ixx,ii) =tmp;
            
            for jj = 1:nGr
                BSr(ii,jj) = mean(double((actual_Yr(:)==jj) - double(predicted_Y(:)==jj)).^2);
            end
            
        end
        BSt =[];
        for jj = 1:nGr
            BSt(jj) =  mean(double((actual_Y(:)==jj) - double(predicted_Y(:)==jj)).^2);
        end
        
        RSSt = 1-repmat(BSt,nBoot,1)./BSr;
        predictedR = [predictedR;nanmean(predictedRt)];
        RSS = [RSS;nanmean(RSSt)];
        BS = [BS;BSt];
        p_sz = [p_sz; nanmean(repmat(BSt,nBoot,1)<BSr)];
        tmp1 = tmp1./nansum(tmp1,2);
        
        
        Cr = nanmean(Cr,3);
        
        IoC = cat(3,IoC,(tmp1-Cr)./Cr);
        FA_rate = sum( (actual_Y==1 & predicted_Y==4))/(sum(actual_Y==1));
        TiW = [TiW;FA_rate];
        
        
        predicted = [predicted;any(predicted_Y(actual_Y==4)==4)];
        
        
        if i ==1
            seizure_ITI = [seizure_ITI;nan];
        else
            seizure_ITI = [seizure_ITI;ix(i)-ix(i-1)];
        end
        
        
        alarmClass = 4;
        refracSec  = 5;     % gaps <= this many seconds are merged into one alarm
        
        % Restrict to true class-1 (interictal) time BEFORE block-merging. isAlarm
        % over the full interval would include the near-onset stretch where a
        % class-4 call is correct, and mixing that into a "false alarm" block count
        % divided by interictal-only seconds double-counts true positives as FAs.
        % This also means a block spanning a class-1 -> class-2 true-label boundary
        % is truncated at the boundary, which is intentional: only the interictal
        % portion counts toward the false-alarm rate.
        interictalMask = actual_Y < 3;
        isAlarm = (predicted_Y == alarmClass) & interictalMask;
        
        % merge alarms separated by gaps <= refracSec, gaps measured in absolute
        % sample index so a run interrupted by non-interictal time is NOT bridged
        % across it
        d = diff([0; isAlarm(:)]);
        onIdx  = find(d == 1);
        offIdx = find(diff([isAlarm(:); 0]) == -1);
        if ~isempty(onIdx)
            keepOn = true(size(onIdx));
            for b = 2:numel(onIdx)
                if onIdx(b) - offIdx(b-1) <= refracSec
                    keepOn(b) = false;   % merge into previous block
                end
            end
            blockOn  = onIdx(keepOn);
            % recompute matching offsets after merge
            keepOff = true(size(offIdx));
            for b = 1:numel(offIdx)-1
                if onIdx(b+1) - offIdx(b) <= refracSec
                    keepOff(b) = false;
                end
            end
            blockOff = offIdx(keepOff);
        else
            blockOn = []; blockOff = [];
        end
        
        nBlockst    = numel(blockOn);
        blockLenMeant = mean(blockOff - blockOn + 1, 'omitnan');
        
        % interictal duration this interval (true class 1), matched to the
        % restriction already applied to isAlarm above
        interictalSec = sum(interictalMask);
        alarmsPerHourt = nBlockst / max(interictalSec,1) * 3600;
        
        % time-in-warning as a fraction, kept for comparison with the old TiW
        TiW_alarmOnlyt = sum(isAlarm) / max(interictalSec,1);
        
        % per-seizure sensitivity: was there >=1 alarm onset in the final 10 s
        % before true onset, i.e. within the true class-4 window. Deliberately
        % NOT restricted to interictalMask -- this is the true positive rate,
        % scored against actual_Y==alarmClass regardless of the class-1 mask
        % used for the false-alarm side above.
        last10 = actual_Y == alarmClass;
        sens_blockt = any(predicted_Y(last10) == alarmClass);
        
        % chance detection rate given the observed block structure of THIS
        % seizure's alarms, for a matched per-seizure null rather than a single
        % pooled number. Treats blocks as a Poisson process at the observed rate
        % with the observed mean length; P(>=1 block overlapping a fixed w-s
        % window) = 1 - exp(-lambda*(L+w)), lambda = blocks/sec, L = mean length,
        % w = horizon length (10 s here).
        w = 10;
        lambda = nBlockst / max(interictalSec,1);
        chance_senst = 1 - exp(-lambda * (blockLenMeant + w));
        alarmsPerHour = [alarmsPerHour;alarmsPerHourt];
        sens_block = [sens_block;sens_blockt];
        TiW_alarmOnly = [TiW_alarmOnly;TiW_alarmOnlyt];
        blockLenMean = [blockLenMean;blockLenMeant];
        chance_sens =[chance_sens;chance_senst];
        
        
        % ================= LATENCY STRATIFICATION OF THE >1 h CLASS =================
        % Splits class 1 (and optionally class 2) by true time-to-seizure so the
        % pooled row in Table 1 can be replaced by something interpretable.
        %
        % Reports per stratum:
        %   rec1      argmax recall (P(predicted==1 | true==1, stratum))
        %   conf1to2  the specific confusion you described: mass going to class 2
        %   conf1oth  mass going to classes 3-6
        %   auc1      one-vs-rest AUC using the class-1 SCORE, not the argmax.
        %             This is the number to report -- it is prevalence-invariant
        %             and independent of the RUSBoost decision threshold.
        
        % ---- before the j / i loops ----
        tedges1 = [3600 7200 14400 28800 Inf];      % 1-2h, 2-4h, 4-8h, >8h
        nS1     = numel(tedges1) - 1;
                         % >= total seizures
        
        rec1     = nan( nS1,1);
        conf1to2 = nan( nS1,1);
        conf1oth = nan( nS1,1);
        auc1     = nan( nS1,1);
        n1       = zeros( nS1,1);
        
        % optional: same treatment for class 2 (100 s - 1 h is itself 1.5 decades)
        tedges2 = [100 316 1000 3600];               % 100-316s, 316-1000s, 1000-3600s
        nS2     = numel(tedges2) - 1;
        rec2 = nan(nS2,1);  auc2 = nan(nS2,1);  n2 = zeros(nS2,1);
        
        % ---- inside the per-seizure loop, REPLACING the existing slice block ----
        % Note this also fixes the confLabel_Y misalignment: the mask is computed
        % once, before actual_Y is overwritten.
        if i == 1
            rng_i = 1:ix(i)+600;
        else
            rng_i = ix(i-1)+600:ix(i)+600;
        end
        
        actual_raw = trueLabel(rng_i);
        pred_raw   = estimateLabel(rng_i, 1);
        conf_raw   = estimateLabel(rng_i, 2:end);
        t2s_raw    = time2seizure(rng_i);
        t2s_raw(t2s_raw>0) = 0;
        
        m           = actual_raw > 0;
        actual_Y    = actual_raw(m);
        predicted_Y = pred_raw(m);
        confLabel_Y = conf_raw(m, :);
        t2s_Y       = abs(t2s_raw(m));       % seconds to the NEXT seizure
        
        % ---- class 1 stratification ----
        s1  = discretize(t2s_Y, tedges1);
        neg = actual_Y ~= 1;                  % shared negative set for AUC
        
        for b = 1:nS1
            pos = (actual_Y == 1) & (s1 == b);
            n1(b) = sum(pos);
            if n1(b) > 0
                rec1(b)     = mean(predicted_Y(pos) == 1);
                conf1to2(b) = mean(predicted_Y(pos) == 2);
                conf1oth(b) = mean(predicted_Y(pos) >  2);
                if any(neg)
                    lab = [true(sum(pos),1); false(sum(neg),1)];
                    sc  = [confLabel_Y(pos,1); confLabel_Y(neg,1)];
                    [~,~,~,auc1(b)] = perfcurve(lab, sc, true);
                end
            end
        end
        
        % ---- optional: class 2 stratification ----
        s2   = discretize(t2s_Y, tedges2);
        neg2 = actual_Y ~= 2;
        for b = 1:nS2
            pos = (actual_Y == 2) & (s2 == b);
            n2(b) = sum(pos);
            if n2(b) > 0
                rec2(i,b) = mean(predicted_Y(pos) == 2);
                if any(neg2)
                    lab = [true(sum(pos),1); false(sum(neg2),1)];
                    sc  = [confLabel_Y(pos,2); confLabel_Y(neg2,2)];
                    [~,~,~,auc2(b)] = perfcurve(lab, sc, true);
                end
            end
        end
        AUC1 = [AUC1;auc1'];
        AUC2 = [AUC2;auc2'];
        N1 = [N1;n1'];
        N2 = [N2;n2'];
        nBlocks = [nBlocks;nBlockst];
    end
    
    
    
    
    
    ID_IoC = [ID_IoC;j*ones(length(tims),1) tims(:)];
    
    all_tim2seiz =[all_tim2seiz;linearize(time2seizure(1:size(estimateLabel,1)))];
    all_pred = [all_pred;estimateLabel];
    all_IDs = [all_IDs;idx*ones(length(estimateLabel),1)];
    
    
    % Merges contiguous class-4 (0-10s) calls into discrete alarms rather than
    % counting raw seconds. refracSec bridges short gaps so one alarm isn't
    % split by a single missed second; adjust to taste.
    
    idx  = idx+1
end




%%


kp = ID_IoC(:,1)~=48 & ID_IoC(:,1)~=18 ;%& (~isnan(seizure_ITI)|seizure_ITI<7200


N = sum(~isnan(p_sz(kp,1)));      % 319
k = sum(p_sz(kp,:) > .95);        % [127 223 190 173 295 303]
pval = 1 - binocdf(k-1, N, 0.05); % upper tail
log10p = log(betainc(0.05, k, N-k+1)) / log(10);   % exact, no underflow
%%
%% ---- repair: zero-alarm seizures have chance 0, not NaN ----
% blockLenMean = mean([]) = NaN when nBlocks == 0, which propagates into
% chance_sens and silently drops those seizures from any nanmean/test.
% Those are the cases with no interictal false alarms at all, i.e. the ones
% where a detection is most informative -- dropping them biases against the
% effect. lambda = 0 there, so chance = 1 - exp(0) = 0.
z = (nBlocks == 0);
blockLenMean(z) = 0;
chance_sens(z)  = 0;

%% ---- assemble the paired vectors ----
ok = kp(:) & ~isnan(sens_block(:)) & ~isnan(chance_sens(:));
obs = logical(sens_block(ok));      % 1 = >=1 alarm in final 10 s
p   = chance_sens(ok);              % matched per-seizure null probability
n   = numel(p);
k   = sum(obs);

fprintf('n seizures        : %d\n', n);
fprintf('observed detected : %d  (%.3f)\n', k, mean(obs));
fprintf('expected (null)   : %.1f  (%.3f)\n', sum(p), mean(p));
fprintf('lift              : %.2f\n', mean(obs)/mean(p));

%% ---- exact Poisson-binomial test (each seizure has its own p_i) ----
pmf = 1;
for ii = 1:n
    pmf = conv(pmf, [1-p(ii), p(ii)]);
end
pval_pb = sum(pmf(k+1:end));        % P(X >= k) under the null
mu  = sum(p);
sd  = sqrt(sum(p.*(1-p)));
zsc = (k - mu)/sd;

fprintf('\nPoisson-binomial  : z = %.1f, p = %.3g\n', zsc, pval_pb);

%% ---- distribution-free cross-check ----
[pval_sr, ~, st] = signrank(double(obs), p, 'tail', 'right');
fprintf('signed-rank       : p = %.3g\n', pval_sr);

%% ---- animal-level test (the one to report; needs session -> rat lookup) ----
% an = rat ID for each seizure, same length/order as sens_block
if exist('an','var')
    ua = unique(an(ok));
    obsA = arrayfun(@(a) mean(obs(an(ok)==a)), ua);
    expA = arrayfun(@(a) mean(p  (an(ok)==a)), ua);
    [~, pval_an] = ttest(obsA, expA, 'Tail', 'right');
    fprintf('\nper-animal (n=%d) : obs %.3f vs exp %.3f, p = %.3g\n', ...
            numel(ua), mean(obsA), mean(expA), pval_an);
    disp([ua(:) obsA(:) expA(:)]);
else
    fprintf('\n[an] not defined -- animal-level test skipped\n');
end
%%
figure
ts = [0 logspace(log10(1),log10(5e4),100)];
imagesc(nanmean(bb,3))

tickVals = [0 10 1e2 1e3 1e4];
xpos = interp1(ts, 1:numel(ts), tickVals);

set(gca, ...
    'XTick', xpos, ...
    'XTickLabel', {'0','10','10^2','10^3','10^4'}, ...
    'TickLabelInterpreter', 'tex', ...
    'XDir', 'reverse');
%%
kp = ID_IoC(:,1)~=48 & ID_IoC(:,1)~=18 ;%& (~isnan(seizure_ITI)|seizure_ITI<7200


figure
plotMeanSEM(1:6,BS(kp,:),'k')

figure
plotMeanSEM(1:6,RSS(kp,:),'k')

figure
plotMeanSEM(1:6,p_sz(kp,:),'k')
%%

IoCt = IoC(:,:,kp);
close all
tmp = nanmean(IoCt*100,3);
figure
imagesc(tmp,[-300 300])
shg
colormap(bluewhitered)
[a,b] = find(tmp>300);
for i = 1:length(a)
    text(b(i)-.25,a(i),num2str(round(10*tmp(a(i),b(i)))/10),...
        'color','w','fontsize',11)
end

set(gca,'ydir','normal')
axis square

colorbar
set(gca,'fontsize',16)

%%
close all
figure
ax  = tight_subplot(6,1);
for i = 1:6
    axes(ax(i))
    plotMeanSEM(1:6,100*squeeze(IoC(:,i,kp))','k')
    hold on
    
    plot([1 6],[0 0])
    plot([5.5 5.5],[10 100],'k')
    set(gca,'xticklabel',[])
    
    set(gca,'fontsize',16)
    xlim([1 6])
    switch i
        case [1,2]
            ylim([-100 100])
        case [3, 4, 6]
            ylim([-200 400])
        case 5
            ylim([-100 10000])
    end
end
%%

close all
figure
imagesc(100*nanmean(C./sum(C,2),3),[0 100])
set(gca,'ydir','normal')
axis square

set(gca,'fontsize',16)
colormap('hot')
colorbar

%%
%sm_ps2pdf(filenameps,filenamepdf,[])
%d = cell2mat(d');

%     % get seizure time series

%     d{idx} = nan(length(tims),220000);
%     for ii = 1:length(tims)
%         d{idx}(ii,:) = LoadBinary([sessions{j,2} '_2.dat'],'frequency',ops.Fs,...
%             'nchannels',ops.nCh_featureFile,'channels',1,'start',tims(ii)-100,'duration',110);
%     end
%     ts = (1:100*ops.Fs)/2000 - 20;
%
%     for ii = 1:length(tims)
%         h=  figure;
%         ax  = tight_subplot(1,2);
%         axes(ax(1))
%         for jj = 1:4
%             tmp = LoadBinary([sessions{j,2} '_' num2str(jj) '.dat'],'frequency',ops.Fs,...
%                 'nchannels',ops.nCh_featureFile,'channels',1,'start',tims(ii)-20,'duration',100);
%             plot(ts,tmp-(jj*300),'k')
%             hold on
%         end
%
%
%         axes(ax(2))
%
%         for jj=1:6
%             plot(ts1,nanconvn(ok1{idx}(ii,:)==jj,k')-(jj*1))
%             hold on
%         end
%         sc = nanmean(ok1{idx}(ii,1900:2000)==4 | ok1{idx}(ii,1900:2000)==3);
%
%         [~,f] = fileparts(sessions{j,1});
%         m = floor(tims(ii)/60);
%         s = floor(mod(tims(ii),60));
%         f = strrep(f,'_',' ');
%         mtit([f '  time: ' num2str(m) ':' num2str(s) ' score '  num2str(round(100*sc)/100)])
%         set(gcf,'position',[1          41        1920        1083])
%         print(h, '-dpsc2',filenameps ,'-append','-bestfit');
%         close all
%
%     end
%
