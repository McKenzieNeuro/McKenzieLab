masterDir = 'R:\McKenzieLab\DANEHippocampalResponse';
fils = getAllExtFiles(masterDir,'mat',1);

%kp_edges = cellfun(@any,regexp(fils,'arena_edges'));
%kp_trans = cellfun(@any,regexp(fils,'transition'));
kp_Novel = (cellfun(@any,regexp(fils,'Novel Env')) | cellfun(@any,regexpi(fils,'LinearTrack\'))) & ~cellfun(@any,regexpi(fils,'Workspaces'));

fils_N = fils(kp_Novel);
[dirs] = fileparts(fils_N);

dirs = unique(dirs);

% get all unique subjects
subj = cellfun(@(a) a(40:44),dirs,'UniformOutput',false);

subj = cellfun(@(a) a(28:32),dirs,'UniformOutput',false);

usubjs = unique(subj);
%%

ts_decay = -.3;
neural =[];
vel =[];
acc = [];
kp_Novel =[];
sesID = [];
ts_trans =[];
subjID =[];

for i = 1:length(dirs)
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        [~,ID] = ismember(subj{i},usubjs);
        
        %load struct with neural/behavioral data
        load('sessiondata.mat')
        if isfield(sessiondata,'contextEntry')
            %get context transitions
            data = sessiondata.contextEntry;
            
            %get times for each context entry
            context_entry = cell2mat(data(:,2));
            
            %sort and get which index goes in which order (b)
            [context_entry,b] = sort(context_entry);
            
            %sort the edge/transition by time
            data = data(b,:);
            
            %define entry times
            epochs_on  = context_entry;
            
            %define exit times
            epochs_off = [context_entry(2:end); sessiondata.behavior.ts_video(end)];
            
            %build matrix of onsets and offsets
            epochs = [epochs_on epochs_off];
            
            %find epochs in home cage
            home = cellfun(@any,regexp(data(:,1),'home'));
            inNovel = ~home;
            
            
            %get times in home cage and in novel context
            ts_video = sessiondata.behavior.ts_video;
            
            %get time from context transition
            
            ts_tran = ts_video;
            
            for k = 1:length(context_entry)-1
                
                kp = ts_video<=context_entry(k+1) & ts_video>context_entry(k);
                ts_tran(kp) = ts_tran(kp) - context_entry(k);
            end
            ts_tran(ts_video>context_entry(end)) = ts_tran(ts_video>context_entry(end)) - context_entry(end);
            
            kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
            kp_homet = cell2mat(kp_homet');
            kp_homet = any(kp_homet,2);
            
            kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
            kp_Novelt = cell2mat(kp_Novelt');
            kp_Novelt = any(kp_Novelt,2);
            
            %get velocity
            velt = sessiondata.behavior.vel;
            acct = sessiondata.behavior.acc;
            distEdgt = sessiondata.behavior.distEdg;
            
            %downsample neural data to video clock
            neuralt = sessiondata.neural.signal_DFoF;
            fs = sessiondata.neural.fs_neural;
            
            
            k = gaussian2Dfilter([1 10*fs],fs);
            neuralt = nanconvn(neuralt,k);
            ts_neural = (1:length(neuralt))/fs;
            
            neuralt = double(interp1(ts_neural,neuralt,ts_video));
            
            nSamples = length(neuralt);
            % concatenate each session
            neural = [neural;neuralt];
            acc = [acc;acct];
            vel = [vel;velt];
            kp_Novel = [kp_Novel;kp_Novelt];
            ts_trans = [ts_trans;ts_tran];
            sesID = [sesID; i *ones(nSamples,1)];
            subjID = [subjID;ID*ones(nSamples,1)];
            i
        end
    end
end
%%


% fit model

%ts_trans1 = ts_trans.^ts_decay;


ix=1;
acc1 = acc;
acc1(acc>75) = 75;
dat = table(neural,vel,acc1,ts_trans,categorical(kp_Novel),categorical(subjID),categorical(sesID), ...
    'VariableNames',{'neural','velocity','acceleration','ts_trans','kp_novel','subjID','sesID'} );

clear vel acc1 acc  sesID sesID ts_trans neural



glme1 = fitglme(dat,'neural ~ 1   + velocity*kp_novel  + (1|sesID:subjID) + (1|subjID)');

glme2 = fitglme(dat,'neural ~ 1   + velocity*kp_novel  + distEdg*kp_novel+  (1|sesID:subjID) + (1|subjID)');

m = compare(glme1,glme2);

LL = glme.LogLikelihood




%%



kp_novel1 = double(table2array(dat(:,5)))==2;
ts_trans1 = table2array(dat(:,4));
vel1 = table2array(dat(:,2));






skewedgaussian = @(amp,alpha,st,offX,offY,x) offY + amp*(1./sqrt((2*pi)).*exp(-(x-offX).^2./(2*st^2))).*normcdf(alpha.*(x-offX));
modelFun =  @(p1,p2,p3,x) p3.* (x./p1).^(p2-1) .* exp(-(x./p1).^p2);
%%

kp_DA =  cellfun(@any,regexp(usubjs(table2array(dat(:,6))),'DA2'));
kp_NE  = ~kp_DA;


dat1 = dat;
neural1 = table2array(dat(:,1));


clear f1 f2 glme resid
for i = 1:2
    
    switch i
        case 1
            kp = kp_DA;
        case 2
            kp = kp_NE;
    end
    ts_trans2 = ts_trans1(kp);
    kp_novel2 = kp_novel1(kp);
    neural2 = neural1(kp);
    
    glme{i} = fitglme(dat(kp,:),'neural ~ 1   + velocity*kp_novel  + (1 |sesID:subjID) + (1 |subjID)');
    ypred = predict( glme{i} );
    
    resid = neural2- ypred;
    
    
    %fit with exponential
    
    
    kp1 = ~isnan(  resid(:,1)) & kp_novel2;
    
    
    f1{i} = fit(ts_trans2(kp1),resid(kp1,1),skewedgaussian, 'Lower',[-100,-10,0,-50,-10],'Upper',[100,10,1000,50,10]);
    
    
    
    
    kp2 = ~isnan(  resid(:,1)) & ~kp_novel2;
    
    f2{i} = fit(ts_trans2(kp2),resid(kp2,1),skewedgaussian, 'Lower',[-100,-10,0,-50,-10],'Upper',[100,10,1000,50,10]);
    
    
    
    
    
    i
end

%%
 close all
for i = 1:2
    
    switch i
        case 1
            kp = kp_DA;
        case 2
            kp = kp_NE;
    end
    ts_trans2 = ts_trans1(kp);
    kp_novel2 = kp_novel1(kp);
    neural2 = neural1(kp);
    ypred = predict( glme{i} );
    resid = neural2- ypred;
    
    
    
    kp = ~isnan(  resid(:,1)) & kp_novel2;
    datBin1  =avghist(ts_trans2(kp),neural2(kp) ,0:1:500);
   
    
    kp = ~isnan(  resid(:,1)) & ~kp_novel2;
    datBin2  =avghist(ts_trans2(kp),neural2(kp) ,0:1:500);
    
   
    
    
    figure
    ax = tight_subplot(3,2);
    ax = reshape(ax,2,3)';
    kp = ~isnan(  resid(:,1)) & kp_novel2;
    pr = f1{i}(ts_trans2(kp));
    axes(ax(1,1))
    plot(0:1:500,datBin1,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),ypred(kp)+pr ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['New Context: data vs full model'])
    
    
    
    axes(ax(2,1))
    plot(0:1:500,datBin1,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),ypred(kp) ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['New Context: data vs linear model'])
    
    axes(ax(3,1))
    plot(0:1:500,datBin1,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),pr ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['New Context: data vs non-linear model'])
    
    kp = ~isnan(  resid(:,1)) & ~kp_novel2;
    
    pr = f2{i}(ts_trans2(kp));
    axes(ax(1,2))
    plot(0:1:500,datBin2,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),ypred(kp)+pr ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['Homecage: data vs full model'])
    
    axes(ax(2,2))
    plot(0:1:500,datBin2,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),ypred(kp) ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['Homecage: data vs linear model'])
    axes(ax(3,2))
    plot(0:1:500,datBin2,'k')
    hold on
    plot(0:1:500,avghist(ts_trans2(kp),pr ,0:1:500),'r')
    ylim([-.25 1.25])
    title(['Homecage: data vs non-linear model'])
    
end
    %%
    
    kp_linear = cellfun(@any,regexp(dirs(table2array(dat(:,7))),'Linear'));
    kp_DA =  cellfun(@any,regexp(usubjs(table2array(dat(:,6))),'DA2'));
    kp_NE  = ~kp_DA;
    %%
    figure
    ax = tight_subplot(2,2);
    for i = 1:4
        axes(ax(i))
        switch i
            case 1
                kp = kp_linear & kp_DA;
                tit = ['DA on linear track'];
            case 2
                kp = kp_linear & kp_NE;
                tit = ['NE on linear track'];
            case 3
                kp = ~kp_linear & kp_DA;
                tit = ['DA in Novel Con'];
            case 4
                kp = ~kp_linear & kp_NE;
                tit = ['NE in Novel Con'];
        end
        kp_novel1 = double(table2array(dat(kp,5)))==2;
        ts_trans1 = table2array(dat(kp,4));
        neural1 = table2array(dat(kp,1));
        neural1 = neural1- ypred(kp);
        tbins = 0:.5:600;
        
        
        plot(tbins,avghist(ts_trans1(kp_novel1),neural1(kp_novel1),tbins),'r')
        hold on
        plot(tbins,avghist(ts_trans1(~kp_novel1),neural1(~kp_novel1),tbins),'k')
        ylim([-1 2])
        title(tit)
    end
    %%
    figure
    
    
    ax = tight_subplot(2,2);
    for i = 1:4
        axes(ax(i))
        switch i
            case 1
                kp = kp_linear & kp_DA;
                tit = ['vel on linear track (DA)'];
            case 2
                kp = kp_linear & kp_NE;
                tit = ['vel on linear track (NE)'];
            case 3
                kp = ~kp_linear & kp_DA;
                tit = ['vel in Novel Con (DA)'];
            case 4
                kp = ~kp_linear & kp_NE;
                tit = ['vel in Novel Con (NE)'];
        end
        kp_novel1 = double(table2array(dat(kp,5)))==2;
        ts_trans1 = table2array(dat(kp,4));
        vel1 = table2array(dat(kp,2));
        tbins = 0:.5:600;
        
        
        plot(tbins,avghist(ts_trans1(kp_novel1),vel1(kp_novel1),tbins),'r')
        hold on
        plot(tbins,avghist(ts_trans1(~kp_novel1),vel1(~kp_novel1),tbins),'k')
        ylim([0 30])
        title(tit)
    end
    
