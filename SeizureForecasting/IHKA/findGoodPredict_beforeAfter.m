


load('R:\IHKA_Scharfman\classification\before_after\classification_2_beforeAfter.mat','sessions')
FeatureFileOutput = 'E:\data\IHKA\features_before_after.mat';
load(FeatureFileOutput,'ops')

filenameps = 'R:\Analysis\McKenzieLab\SeizureForecasting\IHKA\seizSummary_BA1.ps';
filenamepdf = 'R:\Analysis\McKenzieLab\SeizureForecasting\IHKA\seizSummary_BA1.pdf';


%%

group_all = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
dat = cell2mat(dat');
ses = cell2mat(sesID');

lab = [repmat([ops.bins(1:4) ]',4,1) upSample([ ops.bins(5:8)]',4)];
lab = [lab; nan nan];
lab1 = lab(unique(group_all),:);
%%
topDir = 'R:\IHKA_Scharfman\classification\before_after';
fils = getAllExtFiles(topDir,'mat',1);

%%
ix  = 1;C=[];
for i = 2:length(fils)
    tic
    fname = ['R:\IHKA_Scharfman\classification\before_after\classification_' num2str(i) '_beforeAfter.mat'];
    load(fname)
    estimateLabel = predict(rusTree,dat(ses(:,1)==i,:));
    trueLabel = group_all(ses(:,1)==i);
    Ct = confusionmat(trueLabel(trueLabel>0),estimateLabel(trueLabel>0));
    Ctt = nan(17);
    uVal = unique([estimateLabel;trueLabel]);
    Ctt(uVal,uVal) = Ct;
    C(:,:,ix) = Ctt./nansum(Ctt,2);
    ix = ix+1;
    toc
end


%%
figure
[~,b] = sort(lab(:,1));
lab1 = [num2cell(num2str(lab(b(1:end-1),:)),2) ;{'sz'}];
%C1  = C(unique(group_all),unique(group_all),:);
figure
imagesc(nanmean(C(b,b,:),3),[0 .5])
set(gca,'ytick',1:17,'yticklabel',lab1)
set(gca,'xtick',1:17,'xticklabel',lab1,'ydir','normal')
axis square
xlabel('real class')
ylabel('predicted class')

%%
clear dd
for j = 2:16
    
   v= load(['R:\IHKA_Scharfman\prediction\predict_' num2str(j) '_beforeAfter.mat']);
    dd{j} = diff(v.seizure_start);
end
%%
ID =[];C =[];IoC =[];p_sz =[];BS =[];RSS = [];timDiff_post =[];timDiff_pre=[];
idx= 1 ;
clear ok1 d
for j = 2:72
    
    load(['R:\IHKA_Scharfman\prediction\predict_' num2str(j) '_beforeAfter.mat']);
    trueLabel = trueLabel(1:length(estimateLabel));
  
    ts = 1:length(time2seizure);
    Est = nan(17,length(estimateLabel));
    for i = 1:17
        
        Est(i,:) = estimateLabel==i;
        
    end
    [~,ix] =  histc(seizure_start,ts);
    
    ix1 = repmat(ix,1,7201)+repmat(-3600:3600,length(ix),1);
    kp = all(ix1>0 & ix1<size(Est,2),2);
    ix1 = ix1(kp,:);
    ix = ix(kp);
    
    for i = 1:length(ix)
        
        actual_Y = trueLabel(max(ix(i)-3600*6,1):ix(i)+600);
        predicted_Y = estimateLabel(max(ix(i)-3600*6,1):ix(i)+600);
        predicted_Y = predicted_Y(actual_Y>0);
        actual_Y = actual_Y(actual_Y>0);
           ixx  =unique([actual_Y(:);predicted_Y(:)]);
        tmp1 = nan(17,17);
        
        tmp1(ixx,ixx) = confusionmat(actual_Y,predicted_Y);
        
        
        C = cat(3,C,tmp1);
        nBoot = 1000;
        Cr = nan(17,17,nBoot);BSr = nan(nBoot,17);
        for ii = 1:nBoot
            
            % get conf
            actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
            tmp = confusionmat(actual_Yr,predicted_Y);
            
            tmp = tmp./sum(tmp,2);
         
            
            
            Cr(ixx,ixx,ii) =tmp;
            
            for jj = 1:17
                BSr(ii,jj) = mean(double((actual_Yr(:)==jj) - double(predicted_Y(:)==jj)).^2);
            end
            
        end
        BSt =[];
         for jj = 1:17
        BSt(jj) =  mean(double((actual_Y(:)==jj) - double(predicted_Y(:)==jj)).^2);
         end
         
         RSSt = 1-repmat(BSt,1000,1)./BSr;
         
         RSS = [RSS;nanmean(RSSt)];
         BS = [BS;BSt];
         p_sz = [p_sz; nanmean(repmat(BSt,1000,1)<BSr)];
         tmp1 = tmp1./nansum(tmp1,2);
        
         
        Cr = nanmean(Cr,3);
       
        IoC = cat(3,IoC,(tmp1-Cr)./Cr);
       
    end
    
    k  = gaussian2Dfilter([10000 1],10);
    ok1{idx} = zeros(size(ix1));
    
    for i = 1:17
        ok = Est(i,:);
        ok1{idx}= ok1{idx}+ i*ok(ix1);
    end
    
    
    
    %     % get seizure time series
         tims = seizure_start(kp);
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
if length(tims)>1
    timDifft_post = [diff(tims);nan];
        timDifft_pret = [nan;diff(tims)];
else
    timDifft_post = nan(length(tims),1);
    timDifft_pret = nan(length(tims),1);
end

if length(timDifft_post)~=size(ok1{idx},1)
    error('here')
end

    ID = [ID;j*ones(length(tims),1) tims(:)];
    timDiff_post = [timDiff_post;timDifft_post];
     timDiff_pre = [timDiff_pre;timDifft_pret];
    
    idx  = idx+1
end
%sm_ps2pdf(filenameps,filenamepdf,[])
%d = cell2mat(d');
ok1 = cell2mat(ok1');

%%
close all
figure
hold on
y= 100*[nanmean(nanmean([squeeze(IoC(15,16,:)) squeeze(IoC(16,15,:))],2)) nanmean(nanmean([squeeze(IoC(15,15,:)) squeeze(IoC(16,16,:))],2))];
err = 100*[SEM(nanmean([squeeze(IoC(15,16,:)) squeeze(IoC(16,15,:))],2)) SEM(nanmean([squeeze(IoC(15,15,:)) squeeze(IoC(16,16,:))],2))];

bar(y)
errorbar(y,err,'linestyle','none')
figure

y= 100*[nanmean(nanmean([squeeze(IoC(1,5,:)) squeeze(IoC(5,1,:))],2)) nanmean(nanmean([squeeze(IoC(1,1,:)) squeeze(IoC(5,5,:))],2))];
err = 100*[SEM(nanmean([squeeze(IoC(1,5,:)) squeeze(IoC(5,1,:))],2)) SEM(nanmean([squeeze(IoC(1,1,:)) squeeze(IoC(5,5,:))],2))];

bar(y)
hold on
errorbar(y,err,'linestyle','none')


figure


y= 100*[...
   nanmean( nanmean([squeeze(IoC(16,[8],:)) squeeze(IoC(8,16,:))],2)) ...
    nanmean(nanmean([squeeze(IoC(16,[12],:)) squeeze(IoC(12,16,:))],2)) ...
    nanmean(squeeze(IoC(16,16,:)))];
err = 100*[...
   SEM( nanmean([squeeze(IoC(16,[8],:)) squeeze(IoC(8,16,:))],2)) ...
    SEM(nanmean([squeeze(IoC(16,[12],:)) squeeze(IoC(12,16,:))],2)) ...
    SEM(squeeze(IoC(16,16,:)))];


bar(y)
hold on
errorbar(y,err,'linestyle','none')

%%
k  = gaussian2Dfilter([10000 1],3);
for i = 1:size(ok1,1)
    ok3(i,:) = nanconvn(ok1(i,:)==12|ok1(i,:)==16 ,k');
end

figure
[a,b] = max(sum(ok3(:,1900:2000),2),[],2);

[ok2,bb] = sortby(ok1,a);
imagesc(ok2)
ok3 = ok3(bb,:);
%d1 = d(bb,:);

%%
spec = nan(50,22000,size(ok2,1));
for i = 1:size(ok2,1)
    
    sp = abs(awt_freqlist(double(d1(i,:)),2000,logspace(log10(1),log10(250),50))');
    spec(:,:,i) = sp(:,1:10:end);
    i
end


%%
ts = (1:110*ops.Fs/10)/200 - 100;
nopredicted = 1:50;
predicted = 51:size(ok2,1);
close all
freqs = logspace(log10(1),log10(250),50);
spec1 = spec;
spec1(spec>10) = 10;


figure
imagesc(ts,[],nanmean(spec1(:,:,predicted),3),[0 5])
set(gca,'ytick',1:10:50','yticklabel',round(freqs(1:10:end)*10)/10,'ydir','normal')
hold on
plot(ts1,50*nanmean(ok3(predicted,:)),'color',[.7 .7 .7])
xlim([-100 10])
figure
imagesc(ts,[],nanmean(spec1(:,:,nopredicted),3),[0 5])
hold on
plot(ts1,50*nanmean(ok3(nopredicted,:)),'color',[.7 .7 .7])
set(gca,'ytick',1:10:50','yticklabel',round(freqs(1:10:end)*10)/10,'ydir','normal')
xlim([-100 10])

%%
lab = [...

        -Inf          10;...
       -3600          10;...
        -100          10;...
         -10          10;...
        -Inf         100;...
       -3600         100;...
        -100         100;...
         -10         100;...
        -Inf        3600;...
       -3600        3600;...
        -100        3600;...
         -10        3600;...
        -Inf         Inf;...
       -3600         Inf;...
        -100         Inf;...
         -10         Inf;...
         NaN         NaN;...
         
         ];
%%
[~,b] = sort(lab(:,1));

[~,b1] = histc(abs(lab(:,2)),[10 100 3600 inf]);
%b1 = 5-b1;
b1(b1==0) = 5;

kk = timDiff_post<3600;
ok2 = ([flipud(sortby(ok1(~kk,:),timDiff_pre(~kk))) ;((sortby(ok1(kk,:),timDiff_post(kk))))  ]);
ok2a = ok2;
for i = 1:17
    ok2(ok2a ==i) = b1(i);
end

ts1 = (1:size(ok1,2))- 3600;
 col = linspecer(max(ok2(:))-1,'jet');
col = [col ;{ [0 0 0]}];
figure
imagesc(ts1,[],ok2,[.999999999 5.0001])
colormap(cell2mat(col))
colorbar

%%

[~,~,raw] = xlsread('R:\IHKA_Scharfman\IHKA EEG analyses 100321 for Sam.xlsx','SZ analysis');
kp = cellfun(@isstr,raw(:,7));
fname = cell(sum(kp),1);
fname(kp) = cellfun(@(a) strrep(a,' ','_'),raw(kp,7),'uni',0);
% get all seizure classifications
szState =[];szType=[];

TSname1 = 'ure starts';
TSname2 = 'ends';
warning off
N = 0;
for j = 2:21
    seizFil = sessions{j,1};
    TSdata = readtable(seizFil);
    TSdata = table2cell(TSdata);
    seizure_start1 = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname1)),4));
    seizure_start2 = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname1)),4));
    [~,a] = fileparts(sessions{j,2});
  
    ix = find(contains(fname,a));
    ix = ix(ismember(seizure_start2,seizure_start1));
    
    
    load(['R:\IHKA_Scharfman\prediction\predict_' num2str(j) '_beforeAfter.mat']);
    
    ts = 1:length(time2seizure);
    
    [~,ix1] =  histc(seizure_start,ts);
    
    ix1 = repmat(ix1,1,4001)+repmat(-2000:2000,length(ix1),1);
    kp = all(ix1>0 & ix1<length(estimateLabel),2);
    
    
    if length(kp) == length(ix)
        
        ix = ix(kp);
    else
        error('here')
    end
    szState = [szState;raw(ix,20)];
    szType = [szType;raw(ix,23)];
    j
end
%%
szType1= szType(bb);
szState1 =  szState(bb);
 C1 = C./nansum(C,2);
 C1 = C1(:,:,bb);
IoC1 = IoC(:,:,bb);
IoC1 = IoC(:,:,bb);
p_sz1 = p_sz(bb,:);
%%
ts1 = -2000:2000;
col = linspecer(17,'jet')
HYP = contains(szType1,'HYP');
exploration = contains(szState1,'exploration');
sleep = contains(szState1,'Sleep');
immobility = contains(szState1,'immobility');

sleepState = 3*sleep+2*immobility+exploration;
close all
figure
imagesc(ts1,[],ok2(sleepState==3,:),[.999999999 17.0001])
colormap(cell2mat(col))
colorbar

%%
close all
ax  = tight_subplot(3,2);
for i = 1:6
axes(ax(i))
    switch i
        case 1
            kp = HYP & exploration;
        case 2
            kp = HYP & sleep;
        case 3
            kp = HYP & immobility;
        case 4
             kp = ~HYP & exploration;
        case 5
            kp = ~HYP & sleep;
        case 6
             kp = ~HYP & immobility;
    end
    
   
    for ii = 1:17
        for jj = 1:17
            
            [~,pp(ii,jj)] = ttest(squeeze(IoC1(ii,jj,kp)));
        end
    end
    
    d = abs(log10(pp)).*sign(nanmean(IoC1(:,:,kp),3));
    
    imagesc(d,[-30 10])
    colormap(bluewhitered)
    title(num2str(i))
end

%%
close all
figure
ax  = tight_subplot(6,6);
ix=1;
for i = 1:6


    switch i
        case 1
            kp = HYP & exploration;
        case 2
            kp = HYP & sleep;
        case 3
            kp = HYP & immobility;
        case 4
             kp = ~HYP & exploration;
        case 5
            kp = ~HYP & sleep;
        case 6
             kp = ~HYP & immobility;
    end
    t = nanmean(IoC1(:,:,kp),3);
 %  t = t(eye(size(t))==1);
 for j = 1:6
     axes(ax(ix))
     plot(100*t(j,:))
     hold on
     ix = ix+1;
 end
 %ylim([.3 1])
end

%%

close all
figure

for i = 1:3


    switch i
        case 1
            kp =  exploration;
        case 2
            kp =  sleep;
             case 3
            kp = immobility;
      
    end
    t = nanmean(IoC1(:,:,kp),3);
   t = t(eye(size(t))==1);
semilogy(100*t)
hold on
 %ylim([.3 1])
end
legend({'Expl','Sleep','imm'})

%%
figure
ix=1;
for i = 1:2


    switch i
        case 1
            kp =  HYP;
        case 2
            kp =  ~HYP;
          
      
    end
    t = nanmean(IoC1(:,:,kp),3);
    t = t(eye(size(t))==1);
semilogy(100*t)
hold on
end





%%
figure
for ii=1:17
for jj = 1:17
[~,pp(ii,jj)] = ttest2(squeeze(IoC1(ii,jj,HYP)),squeeze(IoC1(ii,jj,~HYP)));
end
end
d = abs(log10(pp)).*sign(nanmean(squeeze(IoC1(:,:,HYP)),3)-nanmean(squeeze(IoC1(:,:,~HYP)),3));
% imagesc(d,[-30 10])
%    colormap(bluewhitered)
 plot(d(eye(size(d))==1))
hold on
 plot(1:17,log10(.05)*ones(17,1),'--')
  plot(1:17,-log10(.05)*ones(17,1),'--')
%%


figure
for ii=1:17
for jj = 1:17
[~,pp(ii,jj)] = ttest2(squeeze(IoC1(ii,jj,sleep)),squeeze(IoC1(ii,jj,~sleep)));
end
end
d = abs(log10(pp)).*sign(nanmean(squeeze(IoC1(:,:,sleep)),3)-nanmean(squeeze(IoC1(:,:,~sleep)),3));
 %imagesc(d,[-30 10])
 %   colormap(bluewhitered)
 plot(d(eye(size(d))==1))
hold on
 plot(1:17,log10(.05)*ones(17,1),'--')
  plot(1:17,-log10(.05)*ones(17,1),'--')
  
  %%
  close all
  for i = 1:17


    switch i
        case 1
            kp = HYP & exploration;
            col = 'k';
            lin = '-';
        case 2
            kp = HYP & sleep;
            col = 'r';
              lin = '-';
        case 3
            kp = HYP & immobility;
               col = 'b';
                 lin = '-';
        case 4
             kp = ~HYP & exploration;
               col = 'k';
                 lin = '--';
        case 5
            kp = ~HYP & sleep;
             col = 'r';
             lin = '--';
        case 6
             kp = ~HYP & immobility;
              col = 'b';
              lin = '--';
    end
    

for ii=1:17
for jj = 1:17
[~,pp(ii,jj)] = ttest(squeeze(IoC1(ii,jj,kp)));
end
end
d = abs(log10(pp)).*sign(nanmean(squeeze(IoC1(:,:,kp)),3));
 %imagesc(d,[-30 10])
 %   colormap(bluewhitered)
 pl(i) = plot(d(eye(size(d))==1),lin,'color',col);
hold on

  end
  
   plot(1:17,log10(.05)*ones(17,1),'k')
  plot(1:17,-log10(.05)*ones(17,1),'k')
  legend(pl,'HYP/Expl','HYP/sleep','HYP/immobility','LVF/Expl','LVF/sleep','LVF/immobility')
  
  %%
  close all
  d=[];
  for i = 1:size(IoC1,3)
  tmp=IoC1(:,:,i);
 tmp = tmp(eye(size(tmp))==1);
d(i,:) = 100*tmp(1:17);
  end
  [a1,bbb]=  sort(nanmean(d(:,3:4),2));
  figure
  plotMeanSEM(1:17,d(sleep& HYP,:),'k')
    plotMeanSEM(1:17,d(~sleep& HYP,:),'r')
    plot(1:17,zeros(1,17),'--')
  ylim([-100 3000])
  %%
  
  
  figure
  plotMeanSEM(1:4,d(~sleep & HYP,:),'k')
    plotMeanSEM(1:4,d(~sleep& ~HYP,:),'r')
      plotMeanSEM(1:4,d(sleep& ~HYP,:),'g')
       plotMeanSEM(1:4,d(sleep& HYP,:),'b')
    plot(1:4,zeros(1,4),'--')
  ylim([-100 3000])
  
  
  %%
  clear d
  close all
  for i = 1:6
      d{i} =  squeeze(IoC1(i,i,:));
   
      
  end
     violin(d)
     
     ylim([-10 300])
     
     %%
     
      close all
  d=[];
  for i = 1:size(IoC1,3)
  tmp=IoC1(:,:,i);
 tmp = tmp(eye(size(tmp))==1);
d(i,:) = 100*tmp;

  end
postictal = nanmean(ok2(:,2000:end)==6,2);
preictal = nanmean(ok2(:,1:2000)==4 | ok2(:,1:2000)==3,2);
  
 pred = nanmean(d(:,[12 16]),2);
  [~,~,~,ix] =   histcn([HYP sleepState],0:1,1:3);
   figure
   pred_state = accumarray(ix+1,postictal,[],@nanmedian,nan);
   pred_state = pred_state(2:end,2:end);
   bar(pred_state)
   
    figure
   pred_state = accumarray(ix+1,pred,[],@nanmedian,nan);
   pred_state = pred_state(2:end,2:end);
   bar(pred_state)
   
   legend('Awake active','Awake inactive','Asleep')
   set(gca,'xtick',1:2,'xticklabel',{'LVF','HYP'})
   
   %%
   close all
   figure
   pred_state = accumarray(ix+1,nanmean(p_sz1(:,3),2)>.99,[],@nanmean,nan);
   pred_state = pred_state(2:end,2:end)
   bar(100*pred_state)
   ylim([0 100])
   figure
   
      pred_state = accumarray(ix+1,nanmean(p_sz1(:,4)>.975 | p_sz1(:,3)>.975 ,2),[],@nanmean,nan);
   pred_state = pred_state(2:end,2:end)
   bar(100*pred_state)
   ylim([0 100])
   %%
   X =[];G=[];
   for i = 1:2
       switch i
           case 1
               kp = HYP==1;
           case 2
               kp = HYP~=1;
       end
       
       for j = 1:3
           kp1 = sleepState==j & kp;
           
           X = [X;pred(kp1)];
           G = [G;((i-1)*3+j)*ones(sum(kp1),1)];
           
       end
   end
   
   
   
   
   boxplot(X,G,'notch','on') 
   
   %%
   
     u = accumarray(ix+1,pred,[],@nanmedian,nan);
     u = u(2:end,2:end)';
     er = accumarray(ix+1,pred,[],@SEM,nan);
     er = er(2:end,2:end)';
     close all
     bar(1:6,u(:),'facecolor','w')
     hold on
      errorbar(1:6,u(:),[],er(:),'linestyle','none','color','k')