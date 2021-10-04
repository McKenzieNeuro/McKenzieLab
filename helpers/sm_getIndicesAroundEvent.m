function [ix,early,late,ts] = sm_getIndicesAroundEvent(evs,pre,post,fs,siz)
%% gets indices around desired events
% INPUTS
%   evs = time of events (s)
%   pre = seconds before
%   post = seconds after
%   fs = sampling rate
%   siz = number of samples in recording
%
% OUTPUT
%   ix = matrix (N = #events, M = time (pre+post)*fs)
%   early = binary vector is any of the rows sample time before 0
%   after = binary vector is any of the rows sample time after 0
%   ts = time stamps of samples around evs

% S. McKenzie, 08/23/2021
%%


evs  = round(fs*evs);
numBack = round(fs*pre);
numForwards = round(fs*post);
totalSample = numBack+numForwards;

ix = repmat(evs(:),1,totalSample) + repmat(-numBack:numForwards-1,length(evs),1);

early = (any(ix<0,2));

late = (any(ix>siz,2));

ts = (-numBack:numForwards-1)/fs;

end