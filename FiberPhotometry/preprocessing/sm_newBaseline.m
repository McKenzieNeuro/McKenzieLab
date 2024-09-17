function sm_newBaseline(dirName,baseEp)


cd(dirName)
if exist([dirName filesep 'sessiondata.mat'])
    
    load([dirName filesep 'sessiondata.mat']);
end


TDTdata = TDTbin2mat(dirName);

% get the photobleached/ isosbestic corrected neural signal (signal_DFoF)
if isfield(TDTdata.streams,'x465A')
    stream = {'x465A','x405A'};
     [signal_DFoF,~,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x405A','baseline',baseEp);
  %use for drug studies where the massive sustained binding alters photobleach correction
  
     % [signal_DFoF,~,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x405A','baseline',baseEp,'photoBleachCorrection','linear','skiptime',baseEp(1),'endTime',baseEp(2));
else
    stream = {'x450D','x500D'};
    %[signal_DFoF,~,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x450D','baseline',baseEp,'photoBleachCorrection','linear','skiptime',baseEp(1),'endTime',baseEp(2));
    [signal_DFoF,~,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x450D','baseline',baseEp);
    
end
sessiondata.neural.signal_DFoF = signal_DFoF;
sessiondata.neural.baseline = baseEp;
save([dirName filesep 'sessiondata.mat'],'sessiondata')
end

