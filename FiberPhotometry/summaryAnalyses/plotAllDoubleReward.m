
function [singles,doubles,other_side,ts_PETH]  = plotAllDoubleReward(dirName)

TDTdata = TDTbin2mat(dirName);

if isfield(TDTdata.epocs,'U12_') &&isfield(TDTdata.epocs,'U11_')
    
    
    
    
    
    
    if isfield(TDTdata.streams,'x465A')
        stream = {'x465A','x405A'};
        [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x405A');
        fs =  TDTdata.streams.x405A.fs;
    else
        stream = {'x450D','x500D'};
        [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x450D');
        fs =  TDTdata.streams.x450D.fs;
    end
    
    %get indices around reward
    [ix2, early, late, ts_PETH] = sm_getIndicesAroundEvent(TDTdata.epocs.U12_.onset,10,10,fs,length(signal_DFoF));
    [ix1, early, late, ts_PETH] = sm_getIndicesAroundEvent(TDTdata.epocs.U11_.onset,10,10,fs,length(signal_DFoF));
    
    
    k = gaussian2Dfilter([100000 1], fs);
    
    % smooth
    ok = nanconvn(signal_DFoF, k');
    %find doubles
    kp_double = find(diff(TDTdata.epocs.U11_.onset)<2);
    nTr = length(TDTdata.epocs.U11_.onset);
    kp_single = setdiff(1:nTr,[kp_double kp_double+1]);
    
    %save for double, single on double side, and single on single side
    singles = ok(ix1(kp_single,:));
    doubles = ok(ix1(kp_double,:));
    other_side = ok(ix2);
else
    
    singles = [];
    doubles =[];
    other_side = [];
    ts_PETH = [];
end

