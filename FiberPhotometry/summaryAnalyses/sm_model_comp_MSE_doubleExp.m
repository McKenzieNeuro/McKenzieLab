function [MSE_full,MSE_red,b_full,pred] = sm_model_comp_MSE_doubleExp(IV,DV,isTime,kp_train,fitRed)
gs = GlobalSearch;

nSamples = size(DV,1);

% do cross val);
kp_test = ~kp_train;

y = DV(kp_train);
y_test = DV(kp_test);

regFactor = 0;
Aeq =[];
beq =[];
A =[];
b0 = [];

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
lb = [-2*ones(1,nLin)];
ub = [2*ones(1,nLin)];





for i = 1:nTime
    d_tmp = avghist(timIV_train(:,i),y,0:20:400);
    eta(i) = (d_tmp(2)-d_tmp(1))/20;
    xo(i) = d_tmp(1);
    xbar(i) = d_tmp(end-1);
    
    
    fitfun_train = [fitfun_train ' + doubleExpFit(timIV_train(:,' num2str(i) ...
        '), b('  num2str(bC+1) '), b('  num2str(bC+2) '),eta(' num2str(i) '), xo(' num2str(i) ...
        '), xbar(' num2str(i) '))'];
    
    fitfun_test = [fitfun_test ' + doubleExpFit(timIV_test(:,' num2str(i) ...
        '), b('  num2str(bC+1) '), b('  num2str(bC+2) '),eta(' num2str(i) '), xo(' num2str(i) ...
        '), xbar(' num2str(i) '))'];
    f = [f ' + b(' num2str(bC+1) ')^2 + b(' num2str(bC+2) ')^2'];
    x0 = [x0 .02 .01];
    lb = [lb .01 .01 ];
    ub = [ub .2 .2];
    bC = bC+2;
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

pred = fitfun_train(b_full);

%[b_full,fg,flg,og] = run(gs,problem);
% lda2 = b_full(6);
% lda1 = b_full(5);
%   C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
%     C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);   
%     xest = C1*exp(-lda1*timIV_train)-C2*exp(-lda2*timIV_train)+xbar;
    
    

%MSE_full = f_test(b_full);%for cross val
figure
yyhat = avghist(timIV_train,pred,1:600);

yy = avghist(timIV_train,y,1:600);

MSE_full = nanmean((yy-yyhat).^2);


% fit reduced model

close all
if fitRed
    for jj = 1:size(IV,2)
        incl = setdiff(1:size(IV,2),jj);
        IV1 = IV(:,incl);
        isTime1 = isTime(incl);
        linearIV_train = IV1(kp_train,~isTime1);
        timIV_train = IV1(kp_train,isTime1);
        
        linearIV_test = IV1(kp_test,~isTime1);
        timIV_test = IV1(kp_test,isTime1);
        
        nLin = size(linearIV_train,2);
        nTime =  size(timIV_train,2);
        fitfun_train = ['@(b) linearIV_train*b(1:' num2str(nLin) ')'''];
        fitfun_test = ['@(b) linearIV_test*b(1:' num2str(nLin) ')'''];
        
        f_train = ['@(b) nanmean((y - fitfun_train(b)).^2 )'];
        f_test = ['@(b) nanmean((y_test - fitfun_test(b)).^2 )'];
        
        f = [f_train ' + ' num2str(regFactor) ' * sum(b(1:' num2str(nLin) ').^2)' ];
        bC = nLin;
        x0 = zeros(1,nLin);
        lb = [-2*ones(1,nLin)];
        ub = [2*ones(1,nLin)];
        
        
        
        for i = 1:nTime
            
            fitfun_train = [fitfun_train ' + doubleExpFit(timIV_train(:,' num2str(i) ...
                '), b('  num2str(bC+1) '), b('  num2str(bC+2) '),eta(' num2str(i) '), xo(' num2str(i) ...
                '), xbar(' num2str(i) '))'];
            
            fitfun_test = [fitfun_test ' + doubleExpFit(timIV_test(:,' num2str(i) ...
                '), b('  num2str(bC+1) '), b('  num2str(bC+2) '),eta(' num2str(i) '), xo(' num2str(i) ...
                '), xbar(' num2str(i) '))'];
            
            %  f = [f ' + b(' num2str(bC+1) ')^2 '];
            x0 = [x0 .01 .001];
            lb = [lb .0001 .0001 ];
            ub = [ub .1 .1];
            bC = bC+2;
        end
        
        
        fitfun_train = eval(fitfun_train);
        fitfun_test = eval(fitfun_test);
        f = eval(f);
        f_train = eval(f_train);
        f_test = eval(f_test);
        
        problem = createOptimProblem('fmincon','objective',f,...
            'x0',x0,'lb',lb,'ub',ub);
        %[b,fg,flg,og] = run(gs,problem);
          b = fmincon(problem);% Objective Function
        
        MSE_red(jj) = f_test(b);
        
    end
else
    MSE_red = nan(size(MSE_full));
end
end

