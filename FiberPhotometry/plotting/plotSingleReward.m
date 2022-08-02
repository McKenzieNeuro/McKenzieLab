function [left,right,tr,ts_PETH]  = plotSingleReward(dirName)


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


    left = ok(ix1);
    right = ok(ix2);
    ts_neural = (1:length(ok))/fs;
    ts = sort([TDTdata.epocs.U11_.onset;TDTdata.epocs.U12_.onset]);
    ts = unique(ts);
    [~,ix] = bestmatch(ts,ts_neural);
    
    tr = nan(length(ix)-1,20000);
    for i = 1:length(ix)-1
    
        tmp = ok(ix(i):ix(i+1));
        tmp = interp1(1:length(tmp),tmp,linspace(1,length(tmp),20000));
        tr(i,:) = tmp;
    end
    
    % get resonse per trial (rewrd to reward)
    
    
else

    left = [];
    ts_PETH = [];
    right = [];
    tr = [];
end
