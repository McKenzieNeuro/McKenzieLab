%load('C:\Users\samckenzie\Downloads\sessionInfo_AAC.mat')
load('R:\AAC\AAC_DataForSam\ArchSessions\AAC\sessionInfo.mat')

ses = 4;

kp = sessionInfo(ses).ispyr;
penalty = mean(sessionInfo(ses).ripBin1,3)/max(mean(sessionInfo(ses).ripBin1,3));
Xtensor = round(sessionInfo(ses).ripBin5(kp,:,:)*.02);

X = sptensor(Xtensor);
%%
% Fit model
opts.K = 20;          % start with 4 templates
opts.maxIter = 300;
opts.lambdaC = 0.05; % optional sparsity on trial weights
opts.verbose = true;
opts.lambda = penalty;
model = poissonCP_trialWeights(X, opts);
%%
X1 = tsne(model.C);

%%


% Xtensor: neurons x time x trials
[nNeurons, nTime, nTrials] = size(Xtensor);

% Sum over time → neurons x trials
Xsum = sum(Xtensor, 2);              % neurons x 1 x trials
Xsum = reshape(Xsum, nNeurons, nTrials);

% Transpose → trials x neurons
Xmat = Xsum.';                       % trials x neurons

[coeff, score, ~] = pca(Xmat, 'NumComponents', 8);
Xmat_PCA = score;   % trials × 8


% Row-normalized version (pattern-only)
rowSums = sum(Xmat, 2);
Xmat_norm = Xmat ./ (rowSums + eps);

rowSums = sum(Xmat_PCA, 2);
Xmat_PCA_norm = Xmat_PCA ./ (rowSums + eps);



Ztensor = model.C;
[coeff, score] = pca(Ztensor);
Ztensor = score(:,1:8);   % manifold space

Ztensor_norm = Ztensor ./ (vecnorm(Ztensor, 2, 2) + eps);
    

%%
optsManifold.idxB = find(sessionInfo(ses).inStim);   % trials to test
optsManifold.idxA = setdiff([1:length(sessionInfo(ses).inStim)]',optsManifold.idxB);   % manifold-defining trials
nB = length(optsManifold.idxB);

optsManifold.k        = 5;
optsManifold.knnType  = 'mutual';    % robust under sparsity
optsManifold.nRefs    = 20;
optsManifold.corrType = 'Spearman';

nShuff = 100;  % number of shuffles for controls

for mo = 1:3
    % ---------------- Select data and options ----------------
    switch mo
        case 1
            rawDat = Xmat;
            
        case 4
            rawDat = Xmat_norm;
           
        case 3
            rawDat = Ztensor;
            optsManifold.transform = 'sqrt';
            optsManifold.patternMetric = 'hcoeffs';
        case 2
            rawDat = Xmat_PCA;
            optsManifold.transform = 'raw';
            optsManifold.patternMetric = 'entropy';
        case 5
            rawDat = Xmat_PCA_norm;
            optsManifold.transform = 'raw';
            optsManifold.patternMetric = 'entropy';
    end

%   result = manifoldAnalysis_ABn(rawDat(optsManifold.idxA,:), ...
 %        rawDat(optsManifold.idxB,:), optsManifold);
    
    % ---------------- Main manifold analysis ----------------
    result(mo).main = manifoldAnalysis_ABn(rawDat(optsManifold.idxA,:), ...
         rawDat(optsManifold.idxB,:), optsManifold);
    % ---------------- Time-shuffle control ----------------
    fracOff_timeShuff = zeros(nB,nShuff);
    fracPattern_timeShuff = zeros(nB,nShuff);
    reconErrB_timeShuff = zeros(nB,nShuff);
    fracOnManifold_timeShuff = zeros(nB,nShuff);
    for s = 1:nShuff
        XtimeShuff = rawDat;

        % Shuffle time bins/features within each trial (row)
        for i = 1:size(XtimeShuff,1)
            XtimeShuff(i,:) = XtimeShuff(i, randperm(size(XtimeShuff,2)));
        end

        % Run manifold analysis
        resTimeShuff = manifoldAnalysis_ABn(rawDat(optsManifold.idxA,:), ...
                                           XtimeShuff(optsManifold.idxB,:), optsManifold);

        fracOff_timeShuff(:,s) = (resTimeShuff.fracOff);
        fracPattern_timeShuff(:,s) = (resTimeShuff.fracPattern);
        reconErrB_timeShuff(:,s) = (resTimeShuff.reconErrB);
        fracOnManifold_timeShuff(:,s) = (resTimeShuff.fracOnManifold);
    end

    result(mo).neuronShuffle.fracOff = fracOff_timeShuff;
    result(mo).neuronShuffle.fracPattern = fracPattern_timeShuff;
    result(mo).neuronShuffle.reconErrB = reconErrB_timeShuff;
result(mo).neuronShuffle.fracOnManifold = fracOnManifold_timeShuff;
end