cd('E:\data\harvard') % now on R drive
clear all
%%
for j = 26:53
    
    %
    if exist(['forecast_summary_' num2str(j) '.mat'])
        load(['forecast_summary_' num2str(j) '.mat'])
        
        trueLabl{j} = [];
        estimateLabel{j} =[];
        sesID{j} = [];
        for i = 1:length(prediction)
            if ~isempty( prediction(i).model)
                prediction1 = prediction(i).model.outpred;
                
                group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),prediction1,num2cell(1:length(prediction1)),'uni',0)');
                prediction1 = cell2mat(prediction1');
                % conf = [conf;cell2mat(prediction(i).model.conf')];
                trueLabl{j} =[trueLabl{j};group];
                estimateLabel{j} = [estimateLabel{j};prediction1];
                sesID{j} = [sesID{j};cell2mat(prediction(i).model.sesID')];
            end
        end
        j
        save('allpred.mat','trueLabl','estimateLabel','sesID','-v7.3')
        clear prediction
    end
    
end
%%
[minSubj,b_min] = min(cellfun(@length,estimateLabel));
ok = cell2mat(cellfun(@(a) a(1:minSubj),estimateLabel,'uni',0));
%ok = cell2mat(estimateLabel);
pred = mode(ok,2);
%m = repmat(pred,1,size(ok,2)) == ok;
%pred(m<.5) = -1;
% tr = trueLabl{b_min};
%  ok1 =  ok(tr==4,:) == 4;
%  [~,b] = max(nanmean(ok1));
%  pr = ok(:,b);
%  sesI = sesID{b_min};
%
[~,b] = unique(sesID{b_min},'rows');

sesI = sesID{b_min}(b,:);
tr = trueLabl{b_min}(b,:);
pr = pred(b,:);
ok = ok(b,:);
%%

clear perf
for i = unique(sesID{1}(:,1))'
    kp = sesI(:,1)==i;
    %  ;
    tr1= tr(kp);
    pr1 = pred(kp);
    perf(i) = nanmean(tr1==pr1);
    %     C = confusionmat(tr1,pr1);
    %     if all(size(C)==[6 6])
    %         C1=C./sum(C,2);
    %         C1(isnan(C1)) = 0;
    %         perf_10s(i,:) = C1(1:4,4);
    %
    %         % get IoC
    %
    %         for j = 1:1000
    %             idx= randsample(1:length(pr1),length(pr1));
    %             pr2 = pr1(idx);
    %             C = confusionmat(tr1,pr2);
    %             C1=C./sum(C,2);
    %             C1(isnan(C1)) = 0;
    %             perf_10sb(j,:) = C1(1:4,4);
    %
    %         end
    %         perf_10s_IoC(i,:) = nanmean((repmat(perf_10s(i,:),1000,1)-perf_10sb)./perf_10sb);
    %     else
    %         perf_10s_IoC(i,:) = nan;
    %         perf_10s(i,:) =nan;
    %     end
    i
end

%%
close all
col = linspecer(10,'jet');
idx=1;
figure
thres = [0 prctile(perf,[10:10:100]) 1];

[~,b]= histc(perf,thres);


for i = 1:10
    
    kp = ismember(sesI(:,1),find(b==i));
    C = confusionmat(tr(kp),pr(kp));
    if all(size(C) == [6 6])
        C1=C./sum(C(:,1:4),2);
        semilogx(-[3*3600 3600 100 10],1-C1(1:4,1),'color',col{i})
        hold on
    end
    i
end


%%
ok = readtable('E:\data\harvard\EEGs_And_Reports_20231024.csv');
gd_ses = sessions(perf>.3,1);
bd_ses = sessions(perf<.3,1);
[~,bidsFolder_gd] = fileparts(fileparts(fileparts(fileparts(gd_ses))));
[~,bidsFolder_bd] = fileparts(fileparts(fileparts(fileparts(bd_ses))));
%%
% loop over subjects

for i = 1:length(gd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_gd{i}))';
    
    str = [lower(ok.reports(idx)); lower(ok.impression(idx))];
    
    hasStroke_G(i) = any(contains(str,'stroke') |  contains(str,'ischem'));
    hasSclerosis_G(i) = any(contains(str,'sclerosis') | contains(str,'lesion'));
    hasTumorStroke_G(i) = any(contains(str,'tumor') | contains(str,'stroke') |  contains(str,'ischem'));
    hasTBI_G(i) =any(contains(str,'tbi')  | contains(str,'traum') & ~contains(str,'no history of trauma') & ~contains(str,'trauma:none'));
    hasMeningitus_G(i) = any(contains(str,'meningitis'));
    hasCongenetial_G(i)= any(contains(str,'congenit') | contains(str,'paulsy'));
    hasGenetic_G(i) = any(contains(str,'genetic') | contains(str,'jme') | contains(str ,'juvenile myoclonic') | contains(str,' gge ') | contains(str,' genetic generalize epilepsy') | contains(str,' ige ') | contains(str,'idiopathatic generalize epilepsy'));% | contains(str,'3hz')  | contains(str,'3 hz'));
    i
    
end

for i = 1:length(bd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_bd{i}))';
    
    str = [lower(ok.reports(idx)); lower(ok.impression(idx))];
    
    hasStroke_B(i) = any(contains(str,'stroke') |  contains(str,'ischem'));
    hasSclerosis_B(i) = any(contains(str,'sclerosis') | contains(str,'lesion'));
    hasTumorStroke_B(i) = any(contains(str,'tumor') | contains(str,'stroke') |  contains(str,'ischem'));
    hasTBI_B(i) =any(contains(str,'tbi')  | contains(str,'traum') & ~contains(str,'no history of trauma') & ~contains(str,'trauma:none'));
    hasMeningitus_B(i) = any(contains(str,'meningitis'));
    hasCongenetial_B(i)= any(contains(str,'congenit') | contains(str,'paulsy'));
    hasGenetic_B(i) = any(contains(str,'genetic') | contains(str,'jme') | contains(str ,'juvenile myoclonic') | contains(str,' gge ') | contains(str,' genetic generalize epilepsy') | contains(str,' ige ') | contains(str,'idiopathatic generalize epilepsy'));% | contains(str,'3hz')  | contains(str,'3 hz'));
    i
    
end


%%

for i = 1:length(bd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_bd{i}))';
    
    
    
    tempSpike_B(i) = any(ok.temporalSpikes(idx)==1);
    occSpike_B(i) =any(ok.occipitalSpikes(idx)==1);
    parSpike_B(i) = any(ok.parietalSpikes(idx)==1);
    
    focalSlowing_B(i) = any(ok.temporalSpikes(idx)==1);
    genSlowing_B(i) =any(ok.occipitalSpikes(idx)==1);
    status_B(i) = any(ok.parietalSpikes(idx)==1);
    K_complexes_B(i) =  any(ok.K_complexes(idx)==1);
    
    i
    
end

for i = 1:length(gd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_gd{i}))';
    
    
    
    tempSpike_G(i) = any(ok.temporalSpikes(idx)==1);
    occSpike_G(i) =any(ok.occipitalSpikes(idx)==1);
    parSpike_G(i) = any(ok.parietalSpikes(idx)==1);
    
    focalSlowing_G(i) = any(ok.temporalSpikes(idx)==1);
    genSlowing_G(i) =any(ok.occipitalSpikes(idx)==1);
    status_G(i) = any(ok.parietalSpikes(idx)==1);
    K_complexes_G(i) =  any(ok.K_complexes(idx)==1);
    
    i
    
end
%%


for i = 1:length(bd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_bd{i}))';
    
    
    
    age_B(i) = mean(ok.AgeAtVisit(idx));
    sex_B(i) = mode(contains(ok.SexDSC(idx),'Female'));
    
    i
    
end

for i = 1:length(gd_ses)
    idx = find(contains(ok.BidsFolder,bidsFolder_gd{i}))';
    
    age_G(i) = mean(ok.AgeAtVisit(idx));
    sex_G(i) = mode(contains(ok.SexDSC(idx),'Female'));
    
    i
    
end



%%
%close all
figure
lab = [{'sclerosis'},{'tumor/stroke'},{'TBI'},{'Infection'},{'genetic'}];
etiol = [...
    
[nanmean(hasSclerosis_G) nanmean(hasSclerosis_B)];
[nanmean(hasTumorStroke_G) nanmean(hasTumorStroke_B)];
[nanmean(hasTBI_G) nanmean(hasTBI_B)];
[nanmean(hasMeningitus_G) nanmean(hasMeningitus_B)];
[nanmean(hasGenetic_G) nanmean(hasGenetic_B)];
];
b = [3  1 4 5  2 ];
%[~,b] = sort((etiol(:,2)-etiol(:,1))./(etiol(:,2)+etiol(:,1)));
barh(etiol(b,:))
set(gca,'ytick',1:5,'yticklabel',lab(b))

%%

[h1,p1] = prop_test([sum(hasSclerosis_G) sum(hasSclerosis_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h2,p2] = prop_test([sum(hasTumorStroke_G) sum(hasTumorStroke_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h3,p3] = prop_test([sum(hasTBI_G) sum(hasTBI_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h4,p4] = prop_test([sum(hasMeningitus_G) sum(hasMeningitus_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h5,p5] = prop_test([sum(hasCongenetial_G) sum(hasCongenetial_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h6,p6] = prop_test([sum(hasGenetic_G) sum(hasGenetic_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
p = [p1 p2 p3 p4 p5 p6];

%%

figure
lab = [{'temp. spike'},{'occ spike'},{'parietal spike'},{'focal slowing'},{'generalized slowing','status','k complex'}];
electro = [...
    
[nanmean(tempSpike_G) nanmean(tempSpike_B)];
[nanmean(occSpike_G) nanmean(occSpike_B)];
[nanmean(parSpike_G) nanmean(parSpike_B)];
[nanmean(focalSlowing_G) nanmean(focalSlowing_B)];
[nanmean(genSlowing_G) nanmean(genSlowing_B)];
[nanmean(status_G) nanmean(status_B)];
[nanmean(K_complexes_G) nanmean(K_complexes_B)];

];

[~,b] = sort((electro(:,2)-electro(:,1))./(electro(:,2)+electro(:,1)));
barh(electro(b,:))
set(gca,'ytick',1:7,'yticklabel',lab(b))

%%

[h1,p1] = prop_test([sum(tempSpike_G) sum(tempSpike_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h2,p2] = prop_test([sum(occSpike_G) sum(occSpike_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h3,p3] = prop_test([sum(parSpike_G) sum(parSpike_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h4,p4] = prop_test([sum(focalSlowing_G) sum(focalSlowing_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h5,p5] = prop_test([sum(genSlowing_G) sum(genSlowing_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h6,p6] = prop_test([sum(status_G) sum(status_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);
[h7,p7] = prop_test([sum(K_complexes_G) sum(K_complexes_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);

p = [p1 p2 p3 p4 p5 p6 p7];
%%
figure
bar(1:2,[nanmean(age_G) nanmean(age_B)],'facecolor','w')
hold on

errorbar(1:2,[nanmean(age_G) nanmean(age_B)],[SEM(age_G') SEM(age_B')],'linestyle','none')

%%
figure
bar(1:2,[nanmean(sex_G) nanmean(sex_B)])
[h8,p8] = prop_test([sum(sex_G) sum(sex_B)] , [length(hasSclerosis_G) length(hasSclerosis_B)], false);

%%

perf1 = perf(perf>0);
GMModel = fitgmdist(perf1',2);
gmPDF = @(x) arrayfun(@(x0) pdf(GMModel,[x0 ]),x);
plot(0:.01:1,gmPDF(0:.01:1)/100)
%%
figure
kp = ismember(sesI(:,1),find(perf>.30));

C = confusionmat(tr(kp),pr(kp));
C1=C./sum(C(:,1:4),2);

kp = ismember(sesI(:,1),find(perf<.30));

C = confusionmat(tr(kp),pr(kp));
C2 = C./sum(C(:,1:4),2);
semilogx(-[5*3600 3600 100 10],C1(1:4,4),'k')
hold on
semilogx(-[5*3600 3600 100 10],C2(1:4,4),'r')
