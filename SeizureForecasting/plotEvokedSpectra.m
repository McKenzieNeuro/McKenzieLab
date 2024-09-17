
function fr = sm_getEventTriggeredSpectra(fname_neural,fname_event,varargin)


%load events
ev = LoadEvents('pulseTimes.evt.sti');

%get first event in series
kp = diff([0;ev.time])>1;
time = ev.time(kp);
%%

freqs = logspace(log10(1),log10(300),100); % frequencies to calc wavelet spectra


pre_event_dur = 10; % time before event
post_event_dur = 10; % time after
dur = pre_event_dur+post_event_dur;


ch = 5; % which channel to load
nCh = 8; % number of channels in binary file
Fs= 20000; % sampling rate
ds_ratio = 100;

fr = nan(length(freqs),(dur*Fs)/ds_ratio,length(time)); %freq x time step x #events


%loop over events
for i = 1:length(time)
    %load channel in window of interest
    
    d = LoadBinary(fname_neural,'frequency',Fs, 'nchannels',nCh,'channels',ch,'start',time(i)-pre_event_dur,'duration',dur);
    
    %get the spectrogram
    fr1 = awt_freqlist(double(d),Fs,freqs);
    
    %down sample because we don't need 20k rez
    fr(:,:,i) = abs(fr1(1:ds_ratio:end,:))';
    
end
end
