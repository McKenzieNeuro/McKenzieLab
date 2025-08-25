% "initialize" variables that will load data across sessions in a loop
%F1-Scoring, sensitivity, specificity, 
%area under the receiver‐operating characteristic curve and time in false warning
ID_IoC =[];C =[];IoC =[];p_sz =[];BS =[];RSS = [];
TiW = []; TiW_r = [];seizure_ITI = [];
idx= 1 ;
%%

%loop over all predictions where each 1-s time window has been categorized
%by which time window it is from

nGr = 6;
clear ok1 d
for j = 1:72

    %load predictions for each 24-hr period
    load(['R:\IHKA_Scharfman\prediction\predict_' num2str(j) '.mat']);
    
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
    kp = all(ix1>0 & ix1<size(Est,2),2);
    ix1 = ix1(kp,:);
    ix = ix(kp);
    tims = seizure_start(kp);
    for i = 1:length(ix)
        if i ==1
        actual_Y = trueLabel(1:ix(i)+600);
        predicted_Y = estimateLabel(1:ix(i)+600);
        else
             actual_Y = trueLabel(ix(i-1)+600:ix(i)+600);
        predicted_Y = estimateLabel(ix(i-1)+600:ix(i)+600);
      
        end
        predicted_Y = predicted_Y(actual_Y>0);
        actual_Y = actual_Y(actual_Y>0);
        ixx  =unique([actual_Y(:);predicted_Y(:)]);
        tmp1 = nan(nGr,nGr);

        tmp1(ixx,ixx) = confusionmat(actual_Y,predicted_Y);

        % we have an observed confusion matrix
        C = cat(3,C,tmp1);
        nBoot = 1000;
        Cr = nan(nGr,nGr,nBoot);BSr = nan(nBoot,nGr);
        for ii = 1:nBoot

            % get conf
            actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
            tmp = confusionmat(actual_Yr,predicted_Y);

            tmp = tmp./sum(tmp,2);



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

        RSS = [RSS;nanmean(RSSt)];
        BS = [BS;BSt];
        p_sz = [p_sz; nanmean(repmat(BSt,nBoot,1)<BSr)];
        tmp1 = tmp1./nansum(tmp1,2);


        Cr = nanmean(Cr,3);

        IoC = cat(3,IoC,(tmp1-Cr)./Cr);
        FA_rate = sum((predicted_Y==3 & actual_Y==1) | (actual_Y==1 & predicted_Y==4))/(sum(actual_Y==1));
        TiW = [TiW;FA_rate];
        
        if i ==1
            seizure_ITI = [seizure_ITI;nan];
        else
            seizure_ITI = [seizure_ITI;ix(i)-ix(i-1)];
        end
    end

    k  = gaussian2Dfilter([10000 1],10);
    ok1{idx} = zeros(size(ix1));

    for i = 1:nGr
        ok = Est(i,:);
        ok1{idx}= ok1{idx}+ i*ok(ix1);
    end


  

    ID_IoC = [ID_IoC;j*ones(length(tims),1) tims(:)];
  

    idx  = idx+1
end

ok1 = cell2mat(ok1');

kp = ID_IoC(:,1)~=48 | ID_IoC(:,1)~=18 ;%& (~isnan(seizure_ITI)|seizure_ITI<7200
%%

figure
plotMeanSEM(1:6,BS(kp,:),'k')
figure
plotMeanSEM(1:6,RSS(kp,:),'k')

figure
plotMeanSEM(1:6,p_sz(kp,:),'k')
%%
[a,b] = sort(sc)
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
