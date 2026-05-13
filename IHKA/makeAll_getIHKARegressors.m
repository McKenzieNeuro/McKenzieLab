
files = getAllExtFiles('R:\DGregg\NeuralData','mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);


%%

% get model per subject and mean CV performance

subjs =  [...
    {'EDS 1.0'} ; ...
    {'EDS 1.3'}; ...
    {'EDS 2.0'}; ...
    {'EDS 2.1'}; ...
    {'EDS 2.2'}; ...
    {'EDS 2.3'}; ...
    {'LTD 10.0'}; ...
    {'EDS 3.0'}; ...
    {'EDS 3.2'}; ...
    {'EDS 4.0'}; ...
    {'EDS 4.1'}; ...
    {'EDS 4.2'}; ...
    {'EDS 5.1'}; ...
    {'EDS 1.1'} ; ...
    {'PCP 4.0'} ; ...
    ];
idx = 1;

clear output
for i = 1:length(subjs)
    kp  = contains(files,subjs{i});
    
    tmp = files(kp);
    
    kp1 = false(size(tmp));
    for j = 1:length(tmp)
        v = load(tmp{j});
        sessiondata = v.sessiondata;
        if isfield(sessiondata,'sleepData')
            kp1(j) = true;
        end
        
    end
    
    tmp = tmp(kp1);
    
    if ~isempty(tmp)
        output{idx} = sm_getIHKARegressors2(tmp);
        idx = idx+1
    end
    
    
end

%%
subjt = cellfun(@(a) a.subj,output,'uni',0);
subj = [];
for i = 1:length(subjt)
    subj = [subj repmat(subjt(i),1,length(output{i}.MSE_base))];
end
    
MSE_sleep = cell2mat(cellfun(@(a) (a.MSE_Sleep-a.MSE_base),output,'uni',0));
MSE_Vel = cell2mat(cellfun(@(a) (a.MSE_Vel-a.MSE_base),output,'uni',0));
MSE_TD = cell2mat(cellfun(@(a) (a.MSE_TD-a.MSE_base),output,'uni',0));
MSE_full= cell2mat(cellfun(@(a) (a.MSE_full-a.MSE_base),output,'uni',0));
MSE_base= cell2mat(cellfun(@(a) (a.MSE_base),output,'uni',0));



  
MSE_sleep_stim = cell2mat(cellfun(@(a) (a.MSE_Sleep_stim-a.MSE_base_stim),output,'uni',0));
MSE_Vel_stim = cell2mat(cellfun(@(a) (a.MSE_Vel_stim-a.MSE_base_stim),output,'uni',0));
MSE_TD_stim = cell2mat(cellfun(@(a) (a.MSE_TD_stim-a.MSE_base_stim),output,'uni',0));
MSE_full_stim= cell2mat(cellfun(@(a) (a.MSE_full_stim-a.MSE_base_stim),output,'uni',0));
MSE_base_stim = cell2mat(cellfun(@(a) (a.MSE_base_stim),output,'uni',0));


%%

kp_subj = ~ismember(subj,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});


figure
bar([nanmean(MSE_full(kp_subj) ) ...
    nanmean(MSE_sleep(kp_subj) ) ...
    nanmean(MSE_TD(kp_subj)) ...
    nanmean(MSE_Vel(kp_subj))])
%%


bd = [nanmean(MSE_full_stim(kp_subj) ) ...
    nanmean(MSE_sleep_stim(kp_subj) ) ...
    nanmean(MSE_TD_stim(kp_subj)) ...
    nanmean(MSE_Vel_stim(kp_subj))];

bs = [SEM(MSE_full_stim(kp_subj)' ) ...
    SEM(MSE_sleep_stim(kp_subj)' ) ...
    SEM(MSE_TD_stim(kp_subj)') ...
    SEM(MSE_Vel_stim(kp_subj)')];
figure
bar(bd,'facecolor','w')
hold on
errorbar(1:4,bd,bs,'linestyle','none','color','k')
set(gca,'xticklabel',{'Full','Sleep','T/D','Vel'})
xlabel('Model')
ylabel('LL ratio (model vs constant)')
%%

vars =  [...
   
    {'SWS'},...
    {'REM'},...
    ];

figure
kp = ~contains(subjt,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
ok  =cell2mat(cellfun(@(a) a.FE_full,output,'UniformOutput',false)');
%[h,p ] =ttest(ok(kp,4))

bd = 100*nanmean(exp(ok(:,3:4))-1);
bs = 100*SEM(exp(ok(:,3:4))-1);
bar(bd,'facecolor','w')
hold on
errorbar(1:2,bd,bs,'linestyle','none','color','k')
set(gca,'xticklabel',vars)
xlabel('Model')
ylabel('Percent in IED rate (from wake)')


%%
hold on
bardat = [ ...
    MSE_full(kp_subj)  ...
    MSE_sleep(kp_subj)  ...
    MSE_TD(kp_subj) ...
    MSE_Vel(kp_subj)];

gr = upSample(1:4,84);

%boxplot(bardat,gr,'notch','on','whisker',inf)
%BoxPlotFormat(gca,'all','Median','Color','k','Whiskers','Color','w')
%%

[~,~,subjt] = unique(subj(kp_subj));
MSE_fullt = MSE_full_stim(kp_subj)';
MSE_sleept = MSE_sleep_stim(kp_subj)';
MSE_TDt = MSE_TD_stim(kp_subj)';
MSE_Velt = MSE_Vel_stim(kp_subj)';
tb = table(MSE_fullt,MSE_sleept,MSE_TDt,MSE_Velt,subjt);
lme = fitlme(tb, 'MSE_sleept ~ 1 + (1|subjt)');



%%
idx=1;
%clear mdl stim
clear stim FE output
for i = 1:length(files)
    v = load(files{i});
    sessiondata = v.sessiondata;
    if isfield(sessiondata,'sleepData')
        
        
        try
            output{idx} = sm_getIHKARegressors(sessiondata);
            
            if sessiondata.StimAmp>100 & contains(lower(sessiondata.StimType),'mono')
                stim(idx) = 1;
            elseif contains(lower(sessiondata.StimType),'none')
                stim(idx) = 0;
            else
                stim(idx) = -1;
            end
            
            
            idx = idx+1;
        end
    end
end
%%
clear badfils
idx = 1;
for i = 1:length(files)
    v = load(files{i});
    sessiondata = v.sessiondata;
    if ~isfield(sessiondata,'sleepData')
        badfils{idx} =      files{i};
        idx = idx+1
    end
end


%%


% get

%%
FE =[];
vars =  [...
    {'TD'},...
    {'vel1'},...
    {'ep_2'},...
    {'ep_3'},...
    {'sleepState1_2'},...
    {'sleepState1_3'},...
    ];

FE = nan(length(mdl),6);
for i = 1:length(mdl)
    [beta,b,c] =   fixedEffects(mdl{i});
    [~,b] = ismember(b.Name,vars);
    
    FE(i,b(b>0)) = beta(b>0);
    
end