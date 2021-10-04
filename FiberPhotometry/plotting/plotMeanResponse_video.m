
%%
%load data
data = TDTbin2mat(pwd);
ts_video = data.epocs.Cam1.onset;
%%

%defin time steps and resample at constant sampling rate

%ts = (1:length(data.streams.data.data)/10)/(data.streams.Wav1.fs/10);

%et 465 sampling rate
fs = data.streams.x465A.fs;

skip = 2e3; %skip file beginning where signal seems to settle
ents = length(data.streams.x465A.data);

%detrend 465
ts_data = (1:length(data.streams.x465A.data))/fs;
a = polyfit(ts_data(skip:ents),data.streams.x465A.data(skip:ents),2);
data_hat465 = polyval(a,ts_data);

%baseline substract 465
data_465= data.streams.x465A.data - data_hat465;

%detrend 405
a = polyfit(ts_data(skip:ents),data.streams.x405A.data(skip:ents),2);
data_hat405 = polyval(a,ts_data);

%baseline substract 405
data_405= data.streams.x405A.data - data_hat405;


% deltaF/F before rescaling?

%data_405 = (data_405-mean(data_405(skip:end)))./mean(data_405(skip:end));

%data_465 = (data_465-mean(data_465(skip:end)))./mean(data_465(skip:end));


%rescale 405
a = polyfit(data_405(skip:ents),data_465(skip:ents),2);

data_405_rescaled = polyval(a,data_405);

%subtract rescaled isobestic
data_NE = data_465-data_405_rescaled;
%%


%calculate deltaF/F after substracting
%data_NE = (data_NE-mean(data_NE(skip:end)))./mean(data_NE(skip:end));

%plot raw data
figure
plot(ts_data,data.streams.x465A.data)

hold on

plot(ts_data,data.streams.x405A.data)

%plot baseline drift
xlim([50 150])
figure
plot(ts_data,data_hat465)

hold on

plot(ts_data,data_hat405)
%plot baseline drift subtracted
xlim([50 150])
figure
plot(ts_data,data_465)

hold on

plot(ts_data,data_405)

%plot with isosbestic rescaled
xlim([50 150])
figure
plot(ts_data,data_465)

hold on

plot(ts_data,data_405_rescaled)

xlim([50 150])

%plot with isosbestic subtracted

figure
plot(ts_data,data_465 - data_405_rescaled)
xlim([50 150])
figure
%plot deltaF/F
plot(ts_data,data_NE)
xlim([50 150])

%%
figure
%pot raw data
plot(ts_data,data.streams.x465A.data)
hold on 
plot (ts_data,data.streams.x405A.data)
%%

%close all
[~,b ] = histc(ts_video(presses),ts_data);

figure


ix = repmat(b(:),1,100001) + repmat(-50000:50000,length(b),1);
d = data_NE(ix); % sample with isosbestic substracted
%d = data_465(ix); % sample from 465 with photbleach subtracted
%d = data.streams.x405A.data(ix); % sample from raw 465

d = nanzscore(d,[],2);


%plot((-5000:5000)/fs,nanmedian(d))

hold on
plot((-50000:50000)/fs,nanmean(d),'k')
%%
figure
imagesc((-50000:50000)/fs,[],d)
%%
figure
ax = tight_subplot(4,5)
ok = nan( 100  ,100001,20);
for i = 1:40
%     axes(ax(i))
% f =logspace(log10(1),log10(10),100);
 wavelet_spectra = awt_freqlist(d(i,:),fs,f);
% 
% imagesc((-50000:50000)/fs,[],abs(wavelet_spectra)')
% set(gca,'ydir','normal','ytick',1:10:100,'yticklabel',round(f(1:10:100)*100)/100)
% xlim([- 10 10])


ok(:,:,i) = abs(wavelet_spectra)';
end

%%

figure
imagesc((-50000:50000)/fs,[],nanzscore(nanmean(ok,3),[],2))
%imagesc((-50000:50000)/fs,[],nanmean(ok,3))
xlim([-5 5])
set(gca,'ydir','normal','ytick',1:10:100,'yticklabel',round(f(1:10:100)*100)/100)
xlim([-10 10])



%%
%get raw data spectrogram
f =logspace(log10(1),log10(10000),100);
d = awt_freqlist(data.streams.Wav1.data(1:10:end),data.streams.Wav1.fs/10,f);
