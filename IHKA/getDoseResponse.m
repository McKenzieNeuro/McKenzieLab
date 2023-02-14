
dirN = 'R:\DGregg\NeuralData\EDS\OL';
fils = getAllExtFiles(dirN,'mat',1);
kp = contains(fils,'recInfo') & contains(fils,'recInfo');

fils = fils(kp);
%%



clear recInfo
[a,b] = fileparts(fils{10});
cd(a)
load('recInfo.mat')
recInfo{3,7} = {'bi','mono'}


save('recInfo.mat','recInfo','-v7.3')


%%

uA = []; mono =[];TD_pre =[];TD_stim =[];subj = [];TD_post =[];
for i = 1:length(fils)
    load(fils{i})
    %get RHS
    RHS = contains(recInfo(:,2),'RHS');
    tmp = recInfo(RHS,:);
    
    uA = [uA tmp{1,4}];
    mono = [mono contains(tmp{1,7},'mono')];
    TD_pre = [TD_pre nanmean(tmp{1,6}{6,2}(1:tmp{1,5}(1))) nanmean(tmp{1,6}{14,2}(1:tmp{1,5}(1)))];
    TD_stim = [TD_stim nanmean(tmp{1,6}{6,2}(tmp{1,5}(1):tmp{1,5}(2))) nanmean(tmp{1,6}{14,2}(tmp{1,5}(1):tmp{1,5}(2)))];
     TD_post = [TD_post nanmean(tmp{1,6}{6,2}(tmp{1,5}(2):end)) nanmean(tmp{1,6}{14,2}(tmp{1,5}(2):end))];
    subj = [subj tmp{1,3}];
end

%%

allsubj = unique(subj);
close all
for i = 1:length(allsubj)
    
    kp = mono==1 & contains(subj,allsubj{i});


subplot(2,1,i)
plot(uA(kp),TD_pre(kp),'o')
hold on
plot(uA(kp),TD_stim(kp),'x')
plot(uA(kp),TD_post(kp),'d')

end

mtit('Monopolar')
figure
for i = 1:length(allsubj)
    
    kp = mono==0 & contains(subj,allsubj{i});


subplot(2,1,i)
plot(uA(kp),TD_pre(kp),'o')
hold on
plot(uA(kp),TD_stim(kp),'x')
plot(uA(kp),TD_post(kp),'d')


end
mtit('Bipolar')