function L = SpkTrainLogLikelihood(q,f,dt)

% L = SpkTrainValuation(S,f,r)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     q: binned spike train
%     f: a tsd describing the predicted intensity function during test epoch
%     dt: time step
% 
% OUTPUT:
%     L: log-likelihood

% Adrien Peyrache, 2014 (following Harris, 2004)

% if length(f) ~= length(q)
%     keyboard
%     error('f and q must be the same length')
% end

L = q.*log(f)-f*dt;
L = sum(L(f>0));
%L = gather(L);

if isinf(L) || isnan(L)
   L =  0;
end