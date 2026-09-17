function [MSE_full,MSE_mat,b_full,pred,MSE_full_tim] = sm_model_comp_MSE_forwards(IV,DV,isTime,kp_train,fitMat)
% Univariate + bivariate model-comparison matrix.
%   MSE_full     : test MSE of the full model (all variables)
%   MSE_mat      : nIV x nIV matrix of test MSEs. Column 1 of IV is the
%                  intercept and is included in every model.
%                    MSE_mat(1,1) = intercept-only baseline
%                    MSE_mat(i,i) = intercept + variable i          (i>=2)
%                    MSE_mat(i,j) = intercept + variables i and j    (i,j>=2, i~=j)
%                    row/col 1 (off-diagonal) = NaN (intercept is not a candidate)
%                  Off-diagonal is symmetric: MSE_mat(i,j) == MSE_mat(j,i).
%   b_full       : full-model coefficients
%   pred         : full-model test-set predictions
%   MSE_full_tim : full-model MSE of the binned (avghist) time course

gs = GlobalSearch;

nSamples = size(DV,1);

% do cross val);
kp_test = ~kp_train;

y = DV(kp_train);
y_test = DV(kp_test);

regFactor = 0.001;
Aeq =[];
beq =[];
A =[];
b0 = [];
bound = 10;


%parse IVs
linearIV_train = IV(kp_train,~isTime);
timIV_train = IV(kp_train,isTime);
linearIV_test = IV(kp_test,~isTime);
timIV_test = IV(kp_test,isTime);

nLin = size(linearIV_train,2);
nTime =  size(timIV_train,2);
fitfun_train = ['@(b) linearIV_train*b(1:' num2str(nLin) ')'''];
fitfun_test = ['@(b) linearIV_test*b(1:' num2str(nLin) ')'''];

f_train = ['@(b) nanmean((y - fitfun_train(b)).^2 )'];
f_test = ['@(b) nanmean((y_test - fitfun_test(b)).^2 )'];
f = [f_train ' + ' num2str(regFactor) ' * sum(b(1:' num2str(nLin) ').^2)' ];
bC = nLin;
x0 = zeros(1,nLin);
lb = [-bound*ones(1,nLin)];
ub = [bound*ones(1,nLin)];

%   eta = (xdata(2)-xdata(1))/dt;
%     xo = xdata(1);
%     C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
%     C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);   
%     xest = C1*exp(-lda1*t)-C2*exp(-lda2*t)+xbar;

for i = 1:nTime
   
    %single
   % fitfun_train = [fitfun_train ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_train(:,' num2str(i) '))'];
   % fitfun_test = [fitfun_test ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_test(:,' num2str(i) '))'];
   
   %     x0 = [x0 0 -.01  ];
%     lb = [lb 0 -.1 ];
%     ub = [ub 20 -.0001];
%     bC = bC+2;
%  f = [f ' + b(' num2str(bC+1) ')^2'];
   
   %double
   
   fitfun_train = [fitfun_train ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_train(:,' num2str(i) ...
       ')) + b(' num2str(bC+3) ')*exp(b(' num2str(bC+4) ')*timIV_train(:,' num2str(i) '))'];
   
    
   fitfun_test = [fitfun_test ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_test(:,' num2str(i) ...
       ')) + b(' num2str(bC+3) ')*exp(b(' num2str(bC+4) ')*timIV_test(:,' num2str(i) '))'];
   f = [f ' + ' num2str(regFactor)  ' * b(' num2str(bC+1) ')^2 +' num2str(regFactor)  '* b(' num2str(bC+3) ')^2'];
  x0 = [x0 0 -.1 0 -.001  ];
    %lb = [lb 0 -.1 -1 -.001 ];
     lb = [lb 0 -.1 -bound -.001 ];
    ub = [ub bound -.001 0 -.0001];
    bC = bC+4;
    
end




beta =[];

% fit full model



fitfun_train = eval(fitfun_train);
fitfun_test = eval(fitfun_test);


f_train = eval(f_train);
f_test = eval(f_test);

f = eval(f);

problem = createOptimProblem('fmincon','objective',f,...
    'x0',x0,'lb',lb,'ub',ub);
b_full = fmincon(problem);

%[b_full,fg,flg,og] = run(gs,problem);
pred =  fitfun_test(b_full);

MSE_full = f_test(b_full);

yyhat = avghist(timIV_test(:,1),pred,1:600);

yy = avghist(timIV_test(:,1),y_test,1:600);

MSE_full_tim = nanmean((yy-yyhat).^2);



% univariate (diagonal) + bivariate (off-diagonal) matrix.
% Column 1 of IV is the intercept and is forced into every model, so it is
% not itself a candidate: the diagonal is intercept+var_i, the off-diagonal
% is intercept+var_i+var_j, row/col 1 are N/A, and MSE_mat(1,1) holds the
% intercept-only baseline.
if fitMat
    nIV     = size(IV,2);
    icept   = 1;                 % intercept column (always included)
    MSE_mat = nan(nIV,nIV);

    % intercept-only baseline
    MSE_mat(icept,icept) = fit_subset(icept,IV,isTime,kp_train,kp_test,y,y_test,regFactor,bound);

    for i = 2:nIV
        % univariate: intercept + variable i
        MSE_mat(i,i) = fit_subset([icept i],IV,isTime,kp_train,kp_test,y,y_test,regFactor,bound);

        % bivariate: intercept + variable i + variable j (upper triangle,
        % mirrored since {i,j} and {j,i} are the same model)
        for j = i+1:nIV
            mse_ij = fit_subset([icept i j],IV,isTime,kp_train,kp_test,y,y_test,regFactor,bound);
            MSE_mat(i,j) = mse_ij;
            MSE_mat(j,i) = mse_ij;
        end
    end
else
    MSE_mat = nan(size(IV,2));
end
end


% -------------------------------------------------------------------------
function [mse_test,mse_train,b] = fit_subset(incl,IV,isTime,kp_train,kp_test,y,y_test,regFactor,bound)
% Fit the linear + double-exponential model on the variable subset `incl`
% and return held-out and training MSE (and the fitted coefficients).

IV1 = IV(:,incl);
isTime1 = isTime(incl);
linearIV_train = IV1(kp_train,~isTime1);
timIV_train    = IV1(kp_train, isTime1);
linearIV_test  = IV1(kp_test, ~isTime1);
timIV_test     = IV1(kp_test,  isTime1);

nLin  = size(linearIV_train,2);
nTime = size(timIV_train,2);

fitfun_train = ['@(b) linearIV_train*b(1:' num2str(nLin) ')'''];
fitfun_test  = ['@(b) linearIV_test*b(1:'  num2str(nLin) ')'''];

f_train = ['@(b) nanmean((y - fitfun_train(b)).^2 )'];
f_test  = ['@(b) nanmean((y_test - fitfun_test(b)).^2 )'];

f  = [f_train ' + ' num2str(regFactor) ' * sum(b(1:' num2str(nLin) ').^2)'];
bC = nLin;
x0 = zeros(1,nLin);
lb = -bound*ones(1,nLin);
ub =  bound*ones(1,nLin);

for i = 1:nTime
    fitfun_train = [fitfun_train ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_train(:,' num2str(i) ...
        ')) + b(' num2str(bC+3) ')*exp(b(' num2str(bC+4) ')*timIV_train(:,' num2str(i) '))'];

    fitfun_test  = [fitfun_test ' + b(' num2str(bC+1) ')*exp(b(' num2str(bC+2) ')*timIV_test(:,' num2str(i) ...
        ')) + b(' num2str(bC+3) ')*exp(b(' num2str(bC+4) ')*timIV_test(:,' num2str(i) '))'];

    f = [f ' + ' num2str(regFactor) ' * b(' num2str(bC+1) ')^2 +' num2str(regFactor) '* b(' num2str(bC+3) ')^2'];

    % exp seed matched to the full model (-.1)
    x0 = [x0 0 -.1 0 -.001];
    lb = [lb 0 -.1 -bound -.001];
    ub = [ub bound -.001 0 -.0001];
    bC = bC+4;
end

fitfun_train = eval(fitfun_train);
fitfun_test  = eval(fitfun_test);
f       = eval(f);
f_train = eval(f_train);
f_test  = eval(f_test);

problem = createOptimProblem('fmincon','objective',f,...
    'x0',x0,'lb',lb,'ub',ub);
b = fmincon(problem);

mse_test  = f_test(b);
mse_train = f_train(b);
end