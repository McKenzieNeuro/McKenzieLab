
load('NE2m3-210819-123801_dataMat.mat')

X = [data_mat(:,3) data_mat(:,4)];
L = @(beta) CostFunction(data_mat(:,2);


options = optimoptions('fmincon','display','off');

weight0 =[];
weight0(1:3:ncat*3) = nbins/2; %offset = 0
weight0(2:3:ncat*3) =  30*(.001/dt); %std = 30ms
weight0([3:3:ncat*3 ncat*3+1:ncat*3+ncont]) =  randn(ncat+ncont,1); %amplidute
%  weight0 = [weight0 rand(1,size(contvar,2))];


problem.lb =[];
problem.lb(1:3:ncat*3) = 1; % pre ms before ripple
problem.lb(2:3:ncat*3) = 10*(.001/dt); % no std
problem.lb([3:3:ncat*3 ncat*3+1:ncat*3+ncont]) = -100*ones(1,ncat+ncont); %random neg amp
% problem.lb = [problem.lb rand(1,size(contvar,2))];

problem.ub =[];
problem.ub(1:3:ncat*3) = nbins; % post  ms before ripple
problem.ub(2:3:ncat*3) = 1000*(.001/dt); % 500ms
problem.ub([3:3:ncat*3 ncat*3+1:ncat*3+ncont]) = 100*ones(1,ncat+ncont); %random neg amp
%  problem.ub = [problem.ub rand(1,size(contvar,2))];


problem.objective = L;
problem.options = options;
problem.solver = 'fmincon';
problem.x0 = weight0;
problem.ObjectiveLimit = 1e-5;


weights = fmincon(problem);