function [home_resp,track_resp,ts_PETH]  = plotTrackResponse(dirName)

if exist('contextTransition.mat')
    v = load('contextTransition.mat');
    ts_home = cell2mat(v.data(cellfun(@any,regexp(v.data(:,1),'home')),2));
    ts_lt = cell2mat(v.data(cellfun(@any,regexp(v.data(:,1),'linear')),2));
else
    error('need to code')
end


if length(ts_home)>1 & ts_home(1)<10
    %assume animal starts in home cage
    ts_home = ts_home(2:end);
end

if ts_lt(1) < 300
    home_resp = [];
    track_resp = [];
    ts_PETH = [];
    return;
    
end


TDTdata = TDTbin2mat(dirName);






if isfield(TDTdata.streams,'x465A')
    stream = {'x465A','x405A'};
    [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x405A');
    fs =  TDTdata.streams.x405A.fs;
else
    stream = {'x450D','x500D'};
    [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x450D');
    fs =  TDTdata.streams.x450D.fs;
end


k = gaussian2Dfilter([100000 1], fs);

% smooth
ok = nanconvn(signal_DFoF, k');

[ix_hc,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts_home,300,300,fs,length(signal_DFoF));
[ix_lt,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts_lt,300,300,fs,length(signal_DFoF));
ok = nanPad(ok,round((ts_home(end)+300) * fs));

home_resp = ok(ix_hc);
track_resp = ok(ix_lt);

end
