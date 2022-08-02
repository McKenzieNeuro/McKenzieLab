function [sharedfit,costval,costval_full] = FitSharedGamma(logISIhist,taubins,varargin)
%FITSHAREDGAMMA Summary of this function goes here
% 
%   INPUT
%       logISIhist      [numtimebins x numcells]  probability density (N/(sum*dbin))
%       taubins         [numtimebins x 1] base e
%   (options)
%       'init_struct'
%       'numAS'         (only needed if no initial guess provided)
%       'AScost_lambda'
%       'AScost_p'
%       'MScost'
%       'MSthresh'
%%
p = inputParser;
addParameter(p,'numAS',3)
addParameter(p,'init_struct',[])
addParameter(p,'AScost_lambda',0.05)
addParameter(p,'AScost_p',1)
addParameter(p,'MScost',10)
addParameter(p,'MSthresh',0.0015)
addParameter(p,'display_results','iter')
addParameter(p,'UseParallel',false)
addParameter(p,'costmean','amean')
addParameter(p,'holdAS',false)

parse(p,varargin{:})
numAS = p.Results.numAS;
init_struct = p.Results.init_struct;

AScost_lambda = p.Results.AScost_lambda;
AScost_p = p.Results.AScost_p;
MScost = p.Results.MScost;
MSthresh = p.Results.MSthresh;
display_results = p.Results.display_results;
UseParallel = p.Results.UseParallel;
costmean = p.Results.costmean;
holdAS = p.Results.holdAS;
%%
numcells = size(logISIhist,2);
%% If there's no initial guess

if isempty(init_struct)
    init_struct.GSlogrates = -0.5.*ones(1,numcells);
    init_struct.GSCVs = 1.5.*ones(1,numcells);
    init_struct.GSweights = 0.5.*ones(1,numcells);

    init_struct.ASlogrates = linspace(1,2.5,numAS);
    init_struct.ASCVs = 0.3.*ones(1,numAS);
    
    init_struct.ASweights  = 0.5.*ones(numcells,numAS)./(numAS);
else
    numAS = length(init_struct.ASlogrates);
end

%%
init = convertGSASparms(init_struct);

%Upper/Lower Bounds
clear lb ub
lb.GSlogrates = -2.*ones(1,numcells);
lb.GSCVs =      zeros(1,numcells);
lb.GSweights =  zeros(1,numcells);
lb.ASlogrates = 0.*ones(1,numAS); %formerly 0.3
lb.ASCVs =      zeros(1,numAS);
lb.ASweights  = zeros(numcells,numAS);

ub.GSlogrates = 2.*ones(1,numcells);
ub.GSCVs =      3.*ones(1,numcells);
ub.GSweights =  ones(1,numcells);
ub.ASlogrates = 3.*ones(1,numAS);
ub.ASCVs =      3.*ones(1,numAS);
ub.ASweights  = ones(numcells,numAS);

if holdAS
    ub.ASlogrates = init_struct.ASlogrates;
    lb.ASlogrates = init_struct.ASlogrates;
    ub.ASCVs = init_struct.ASCVs;
    lb.ASCVs = init_struct.ASCVs;
end

lb = convertGSASparms(lb);
ub = convertGSASparms(ub);

% typ.GSlogrates = -1.*ones(1,numcells);
% typ.GSCVs =      1.*ones(1,numcells);
% typ.GSweights =  0.5.*ones(1,numcells);
% typ.ASlogrates = 1.*ones(1,numAS);
% typ.ASCVs =      0.3.*ones(1,numAS);
% typ.ASweights  = 0.1.*ones(numcells,numAS);
% typ = convertGSASparms(typ);

%Make the constraint matrix for all weights to add to 1
Aeq = zeros(numcells,length(ub));
Aeq_ASonly = zeros(numcells,length(ub));
Beq = ones(numcells,1);
for cc = 1:numcells
    thiscell.GSlogrates = zeros(1,numcells);
    thiscell.GSCVs =      zeros(1,numcells);
    thiscell.GSweights =  zeros(1,numcells);
    thiscell.ASlogrates = zeros(1,numAS);
    thiscell.ASCVs =      zeros(1,numAS);
    thiscell.ASweights  = zeros(numcells,numAS);
    thiscell.ASweights(cc,:) = 1;
    Aeq_ASonly(cc,:) = convertGSASparms(thiscell);
    thiscell.GSweights(cc) = 1;
    Aeq(cc,:) = convertGSASparms(thiscell);
end
Aeq_ASonly(Aeq_ASonly~=1)=0;
Aeq(Aeq~=1)=0;

options = optimoptions('fmincon','Algorithm','sqp' ,'UseParallel',UseParallel,'Display',display_results);%
%try also: 'Algorithm','interior-point''active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e8;
options.MaxIterations = 1500; 
options.StepTolerance = 1e-8;
% options.TypicalX = typ;
%options.FiniteDifferenceStepSize = 1000.*sqrt(eps);

%% Fit all the distributions together

%The Loss Function for each cell
cellloss = @(GSASparm_vect) sum((logISIhist-GSASmodel2(GSASparm_vect,taubins,numcells,numAS)).^2).^0.5;

%Loss function for only refreactory spikes
%Only for high spike density (positive)
smallISI = MSthresh; %<2ms ISIs penalize (make parameter)
sub1msbins = taubins<=log(smallISI); %Which bins are small enough for small-time cost
cellloss_ref = @(GSASparm_vect) sum(...
    (logISIhist(sub1msbins,:)-GSASmodel2(GSASparm_vect,taubins(sub1msbins),numcells,numAS)).^2).^0.5;

switch costmean
    case 'gmean'
        %Total loss function with regularization etc
        costfun = @(GSASparm_vect) 10.^((1./numcells).*sum(log10(... %Geometric mean (so small # bad cells don't dominate)
            cellloss(GSASparm_vect) + ...
            AScost_lambda.*(abs(Aeq_ASonly*GSASparm_vect)').^(AScost_p) + ...; %L1/2 norm on AS weights to promote sparseness
            MScost.*cellloss_ref(GSASparm_vect)))); 

    case 'amean'
        %Total loss function with regularization etc
        costfun = @(GSASparm_vect) (1./numcells).*sum(... %Geometric mean (so small # bad cells don't dominate)
            cellloss(GSASparm_vect) + ...
            AScost_lambda.*(abs(Aeq_ASonly*GSASparm_vect)').^(AScost_p) + ...; %L1/2 norm on AS weights to promote sparseness
            MScost.*cellloss_ref(GSASparm_vect)); 
end

%%
%Fitting
fitparms = fmincon(costfun,init,[],[],Aeq,Beq,lb,ub,[],options);
costval_full = costfun(fitparms); %Get the total loss
costval = cellloss(fitparms); %Get the loss for each cell

%Convert back into structure for output
sharedfit = convertGSASparms(fitparms,numcells,numAS);

%Sort AS modes by mean rate from high to low
[~,ASratesort] = sort(sharedfit.ASlogrates,'descend');
sharedfit.ASlogrates = sharedfit.ASlogrates(ASratesort);
sharedfit.ASCVs = sharedfit.ASCVs(ASratesort);
sharedfit.ASweights = sharedfit.ASweights(:,ASratesort);

end

