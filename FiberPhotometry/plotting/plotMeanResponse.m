
cd E:\data\Donaldson\NE2m3-210518-101716
cd E:\data\Donaldson\NE2h2-210518-105827
cd E:\data\Donaldson\NE2h1-210518-113529
cd R:\DonaldsonT\NE2h1-210518-113529

%%
data = TDTbin2mat(pwd);
%%

f =logspace(log10(1),log10(10000),100);
d = awt_freqlist(data.streams.Wav1.data(1:10:end),data.streams.Wav1.fs/10,f);
ts = (1:length(data.streams.Wav1.data)/10)/(data.streams.Wav1.fs/10);
%%
oddball = 35; %23

expected = 62; %291 
figure
plot(abs(d(:,expected)),'k')
hold on
plot(abs(d(:,oddball))) %oddball

%%





%find oddball
odd_thres = .02;
com_thres = .005;

ts_odd = ts(find(diff(abs(d(:,oddball))>odd_thres)>0));



ts_com = ts(find(diff(abs(d(:,expected))>com_thres)>0));


%%
skip = 1e4;
fs = data.streams.x465A.fs;
ts_data = (1:length(data.streams.x465A.data))/fs;
a = polyfit(ts_data(skip:end),data.streams.x465A.data(skip:end),2);

data_hat = polyval(a,ts_data);

%baseline substract

data_465= data.streams.x465A.data - data_hat;



a = polyfit(ts_data(skip:end),data.streams.x405A.data(skip:end),2);

data_hat = polyval(a,ts_data);

%baseline substract

data_405= data.streams.x405A.data - data_hat;




%without isobestic
%data_465 = (data_465-median(data_465(skip:end)))./median(data_465(skip:end));
%data_405 = (data_405-median(data_405(skip:end)))./median(data_405(skip:end));

a = polyfit(data_405(skip:end),data_465(skip:end),1);
data_405_rescale = polyval(a,data_405);

data_NE = data_465-data_405_rescale;
%data_NE = (data_NE-nanmedian(data_NE(skip:end)))/nanmedian(data_NE(skip:end));
%%
%plot signal


[~,b ] = histc(ts_odd,ts_data);

[~,b1 ] = histc(ts_com,ts_data);


ix = repmat(b(:),1,10001) + repmat(-5000:5000,length(b),1);

ix1 = repmat(b1(:),1,10001) + repmat(-5000:5000,length(b1),1);

d1 = data_NE(ix1);
d = data_NE(ix);
d = nanzscore(d,[],2);

d1 = nanzscore(d1,[],2);
figure
plot((-5000:5000)/fs,nanmean(d))

hold on
plot((-5000:5000)/fs,nanmean(d1),'k')

%%
f =logspace(log10(5),log10(1000),100);
data_NE_freq =awt_freqlist(data_NE,fs,f);
ok = nan(10001,100,900);
for i = 1:size(ix1,1)
   
    ok(:,:,i) = abs(data_NE_freq(ix1(i,:),:));
    
    
    
end

%%

u = nanmedian(ok,3)';
figure
imagesc(nanzscore(nanmedian(ok,3)',[],2))
set(gca,'ydir','normal')
