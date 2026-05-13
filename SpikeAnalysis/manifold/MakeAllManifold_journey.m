clear all
topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\AAC\Acute';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\PV';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\CCK';
plotIt = false;
d = dir(topDir);
d = {d(cell2mat({d.isdir})).name};
dirs = d(3:end);

dirs = cellfun(@(a) [topDir filesep a],dirs,'uni',0);
kp = ~contains(dirs,'figures');
dirs = dirs(kp);
clear p_comp uRate_comp spont_order stim_order
%%
% loop over directories

for se = 1:length(dirs)

    cd(dirs{se})
    fils = getAllExtFiles(pwd,'mat',0);

    kp = contains(fils,'celltypes') |  contains(fils,'ripples') |  contains(fils,'optoStim')  |  contains(fils,'spikes') | contains(fils,'cell_metrics');
    fils = fils(kp);
    %load all data

    %load data
    clear allcelltypes ripples spikes optoStim cell_metrics
    for j = 1:length(fils)
        load(fils{j})
    end

    if exist('allcelltypes')
        ispyr =cellfun(@(a) isstr(a) && contains(a,'pyr'),allcelltypes);

    else
        ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);
    end
    in = InIntervals(ripples.peaks,optoStim.timestamps);

    [binnedPop5,bin_times]=populationMatrix(spikes,0,.1,5,ripples.timestamps(:,1));
    [binnedPop1,bin_times]=populationMatrix(spikes,0,.1,1,ripples.timestamps(:,1));

    sessionInfo(se).name = dirs{se};
    sessionInfo(se).ripple_times = ripples.timestamps;
    sessionInfo(se).stim_times = optoStim.timestamps;
    sessionInfo(se).inStim = in;
    sessionInfo(se).ripBin1 = binnedPop1;
    sessionInfo(se).ripBin5 = binnedPop5;
    sessionInfo(se).ispyr=ispyr;



    kp = sessionInfo(se).ispyr;

    % Transpose → trials x neurons
    Xmat = squeeze(sessionInfo(se).ripBin1(kp,:,:));                      % trials x neurons
    [coeff, score, ~] = pca(Xmat', 'NumComponents', 8);
    Xmat_PCA = score;   % trials × 8

    optsManifold.idxB = find(sessionInfo(se).inStim);   % trials to test
    optsManifold.idxA = setdiff([1:length(sessionInfo(se).inStim)]',optsManifold.idxB);   % manifold-defining trials
    rIDXA = optsManifold.idxA(randsample(length(optsManifold.idxA),length(optsManifold.idxA)));
    nBlock = 10;
    nShuff = 10;  % number of shuffles for controls
    blidx = 1:floor(length(rIDXA)/nBlock):length(rIDXA);
    [n,bl] = histc(1:length(rIDXA),blidx);

    %define k-fold
    for b = 1:nBlock

        manA = rIDXA(bl ~=b);
        testA = rIDXA(bl ==b);
        nB = length(optsManifold.idxB);

        optsManifold.k        = 5;
        optsManifold.knnType  = 'mutual';    % robust under sparsity
        optsManifold.nRefs    = 20;
        optsManifold.corrType = 'Spearman';



        rawDat = Xmat';
        %rawDat = Xmat_PCA;
        result(se,b).mainA = manifoldAnalysis_ABn( rawDat(manA,:),...
            rawDat(testA,:), optsManifold);

        result(se,b).mainB = manifoldAnalysis_ABn( rawDat(manA,:), ...
            rawDat(optsManifold.idxB,:), optsManifold);

        % ---------------- Time-shuffle control ----------------
        fracOff_timeShuff = zeros(nB,nShuff);
        fracPattern_timeShuff = zeros(nB,nShuff);
        reconErrB_timeShuff = zeros(nB,nShuff);
        fracOnManifold_timeShuff = zeros(nB,nShuff);
        fracGain_timeShuff = zeros(nB,nShuff);
        for s = 1:nShuff
            XtimeShuff = rawDat;

            % Shuffle time bins/features within each trial (row)
            for i = 1:size(XtimeShuff,1)
                XtimeShuff(i,:) = XtimeShuff(i, randperm(size(XtimeShuff,2)));
            end

            % Run manifold analysis
            resTimeShuff = manifoldAnalysis_ABn(rawDat(manA,:), ...
                XtimeShuff(optsManifold.idxB,:), optsManifold);

            fracOff_timeShuff(:,s) = (resTimeShuff.fracOff);
            fracPattern_timeShuff(:,s) = (resTimeShuff.fracPattern);
            reconErrB_timeShuff(:,s) = (resTimeShuff.reconErrB);
            fracOnManifold_timeShuff(:,s) = (resTimeShuff.fracOnManifold);
            fracGain_timeShuff(:,s) = (resTimeShuff.fracOnManifold);

        end

        result(se,b).neuronShuffle.fracOff = fracOff_timeShuff;
        result(se,b).neuronShuffle.fracPattern = fracPattern_timeShuff;
        result(se,b).neuronShuffle.reconErrB = reconErrB_timeShuff;
        result(se,b).neuronShuffle.fracOnManifold = fracOnManifold_timeShuff;
        result(se,b).neuronShuffle.fracGain = fracGain_timeShuff;

    end
end
%%

fracOff_randNeuron = mean(arrayfun(@(a) mean(mean(a.neuronShuffle.fracOff)),result),2);
fracOff_mainA = mean(arrayfun(@(a) mean(mean(a.mainA.fracOff)),result),2);
fracOff_mainB = mean(arrayfun(@(a) mean(mean(a.mainB.fracOff)),result),2);

fracGain_randNeuron = mean(arrayfun(@(a) mean(mean(a.neuronShuffle.fracGain)),result),2);
fracGain_mainA = mean(arrayfun(@(a) mean(mean(a.mainA.fracGain)),result),2);
fracGain_mainB = mean(arrayfun(@(a) mean(mean(a.mainB.fracGain)),result),2);


fracPattern_randNeuron = mean(arrayfun(@(a) mean(mean(a.neuronShuffle.fracOnManifold)),result),2);
fracPattern_mainA = mean(arrayfun(@(a) mean(mean(a.mainA.fracOnManifold)),result),2);
fracPattern_mainB = mean(arrayfun(@(a) mean(mean(a.mainB.fracOnManifold)),result),2);

reconErrB_randNeuron = mean(arrayfun(@(a) mean(mean(a.neuronShuffle.reconErrB)),result),2);
reconErrB_mainA = mean(arrayfun(@(a) mean(mean(a.mainA.reconErrB)),result),2);
reconErrB_mainB = mean(arrayfun(@(a) mean(mean(a.mainB.reconErrB)),result),2);

reconErrB_randNeuron = mean(arrayfun(@(a) mean(mean(a.neuronShuffle.reconErrB)),result),2);
reconErrB_mainA = mean(arrayfun(@(a) mean(mean(a.mainA.reconErrB)),result),2);
reconErrB_mainB = mean(arrayfun(@(a) mean(mean(a.mainB.reconErrB)),result),2);


%%
figure
ax = tight_subplot(3,1);
axes(ax(1))
X = mean([fracOff_mainA fracOff_mainB fracOff_randNeuron]);
S = SEM([fracOff_mainA fracOff_mainB fracOff_randNeuron]);
bar(X,'w')
hold on
errorbar(1:3,X,S,'LineStyle','none','color','k')
ylabel('Fraction off manifold')
set(gca,'xtick',1:3,'xticklabel','','fontsize',16)
%axes(ax(2))

%bar(mean([fracGain_mainA fracGain_mainB fracGain_randNeuron]))


axes(ax(2))
X = mean([reconErrB_mainA reconErrB_mainB reconErrB_randNeuron]);
S = SEM([reconErrB_mainA reconErrB_mainB reconErrB_randNeuron]);
bar(X,'w')
hold on
errorbar(1:3,X,S,'LineStyle','none','color','k')
ylabel('Reconstruction error')

set(gca,'xtick',1:3,'xticklabel',{'w/AAC','w/o AAC','shuffle'},'fontsize',16)


axes(ax(3))

bar(mean([fracPattern_mainA fracPattern_mainB fracPattern_randNeuron]))