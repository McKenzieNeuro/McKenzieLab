function [home_resp1,home_resp2,track_respNO,track_respLT1,ts_PETH]  = plotTrackResponseNovel(dirName)
cd(dirName)
if exist('contextTransition.mat')
    v = load('contextTransition.mat');
    ts_home = cell2mat(v.data(cellfun(@any,regexp(v.data(:,1),'home')),2));

    NO = cellfun(@any,regexp(v.data(:,1),'linear_track_NO'));
    lt = cellfun(@any,regexp(v.data(:,1),'linear'));
    lt1 = xor(lt,NO);
    
    ts_NO = cell2mat(v.data(lt,2));
    ts_lt1 = cell2mat(v.data(lt1,2));
else
    error('need to code')
end


if length(ts_home)>1 & ts_home(1)<10
    %assume animal starts in home cage
    ts_home = ts_home(2:end);
end

if ts_lt1(1) < 300
   home_resp = [];
   track_respNO =[];
   track_respLT1 =[];
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

[ix_hc1,early,late,ts_PETH{1}] = sm_getIndicesAroundEvent(ts_home(1),75,30,fs,length(signal_DFoF));
[ix_hc2,early,late,ts_PETH{2}] = sm_getIndicesAroundEvent(ts_home(2),300,300,fs,length(signal_DFoF));
[ix_ltNO,early,late,ts_PETH{3}] = sm_getIndicesAroundEvent(ts_NO,30,300,fs,length(signal_DFoF));
[ix_ts_lt1,early,late,ts_PETH{4}] = sm_getIndicesAroundEvent(ts_lt1,300,75,fs,length(signal_DFoF));


ok = nanPad(ok,round((ts_home(end)+300) * fs));

home_resp1 = ok(ix_hc1);
home_resp2 = ok(ix_hc2);
track_respNO = ok(ix_ltNO);
track_respLT1 = ok(ix_ts_lt1);
end
