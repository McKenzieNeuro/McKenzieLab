fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexpi(fils,'GEFI'));
fils  =fils(kp);



%%
filenameps = 'R:\DANEHippocampalResponse\linearTackSummary.ps';
filenamepdf = 'R:\DANEHippocampalResponse\linearTackSummary.pdf';
evFile = 'trackON.mat';
for i = 1:length(fils)
    
    
    sm_plotGEFI_Summary(fils{i},evFile,filenameps)
    close all
end


ps2pdf('psfile',filenameps,'pdffile',filenamepdf,'gscommand','C:\Program Files (x86)\gs\gs9.54.0\bin\gswin32.exe')



%%
fs = 1017.3;ff=[];
clear track
evFile = 'trackON.mat';
for i = 1:length(fils)
    
    
    k  = gaussian2Dfilter([1 10000],fs);
    [a,b] = fileparts(fils{i});
    cd(a)
    load(fils{i})
    
    
    isNE = regexp(a,'NE2');
    isDA = regexp(a,'DA2');
    
    if isempty(isDA) & isempty(isNE)
        isDA = regexp(a,'DACSD');
    end
    
    
    
    if ~isempty(evFile)
        e = load([a filesep evFile]);
    end
    
    
    
    tr = cell2mat(e.data(cellfun(@any,regexp(e.data(:,1),'track'))& ~cellfun(@any,regexp(e.data(:,1),'track_')),2));
    trN = cell2mat(e.data(cellfun(@any,regexp(e.data(:,1),'track_')),2));
    HC = cell2mat(e.data(cellfun(@any,regexp(e.data(:,1),'home')),2));
    
    [ix1, early, late, ts_PETH] = sm_getIndicesAroundEvent(GEFI.reward1,60,60,GEFI.fs,length(GEFI.corrected));
    [ix2, early, late, ts_PETH] = sm_getIndicesAroundEvent(GEFI.reward2,60,60,GEFI.fs,length(GEFI.corrected));
    [ix3, early, late, ts_PETH] = sm_getIndicesAroundEvent(tr,300,300,GEFI.fs,length(GEFI.corrected));
    
    [ix4, early, late, ts_PETH] = sm_getIndicesAroundEvent(HC,300,300,GEFI.fs,length(GEFI.corrected));
    
    
    ix = [ix1;ix2];
    
    signal_smoothed = nanconvn(GEFI.corrected, k);
    [f,tmp] = plotFFT(signal_smoothed,fs,0);
    ff = [ff;avghist(f,double(tmp),0:.01:5)];
    reward{i} = signal_smoothed(ix);
    track{i} = signal_smoothed(ix3);
    homecage{i} = signal_smoothed(ix4);
    
    if~isempty(trN)
        
        [ix5, early, late, ts_PETH] = sm_getIndicesAroundEvent(trN,300,300,GEFI.fs,length(GEFI.corrected));
        trackN{i} = signal_smoothed(ix5);
    else
        trackN{i} = nan(1,length(ts_PETH));
    end
    
    
    indicator(i) = any(isDA);
    power_signal(i,:) = [GEFI.power_isos GEFI.power_signal];
    i
end

%%

kpDA = find(power_signal(:,2)<20 & indicator(:)==1);
kpNE = find(power_signal(:,2)<10 & indicator(:)==0);
%%
close all
figure
imagesc(ts_PETH,[],(zscore(cell2mat(reward(kpNE)'),[],2)),[-3 3])
figure
imagesc(ts_PETH,[],(zscore(cell2mat(reward(kpDA)'),[],2)),[-3 3])

%%

rewardu = cell2mat(cellfun(@(a) mean(a,1),reward,'UniformOutput',false)');
HCu = cell2mat(cellfun(@(a) mean(a,1),homecage,'UniformOutput',false)');
tracku = cell2mat(cellfun(@(a) mean(a,1),track,'UniformOutput',false)');
trackNu = cell2mat(cellfun(@(a) nanmean(double(a),1),trackN,'UniformOutput',false)');
%%
figure
subplot(1,3,1)
imagesc(ts_PETH,[],trackNu(kpDA,:),[-1 1])

subplot(1,3,2)
imagesc(ts_PETH,[],tracku(kpDA,:),[-1 1])


subplot(1,3,3)
imagesc(ts_PETH,[],HCu(kpDA,:),[-1 1])
%%
figure
subplot(1,3,1)
imagesc(ts_PETH,[],trackNu(kpNE,:),[-1 1])

subplot(1,3,2)
imagesc(ts_PETH,[],tracku(kpNE,:),[-1 1])


subplot(1,3,3)
imagesc(ts_PETH,[],HCu(kpNE,:),[-1 1])


