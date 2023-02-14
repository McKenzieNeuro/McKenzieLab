


d1 = LoadBinary('AC75a-5_DOB_072519_TS_2020-03-23_17_30_04_1.dat','nchannels',41,'channels',1,'frequency',2000,'duration',17000)
d2 = LoadBinary('AC75a-5_DOB_072519_TS_2020-03-23_17_30_04_2.dat','nchannels',41,'channels',1,'frequency',2000,'duration',17000)
d3 = LoadBinary('AC75a-5_DOB_072519_TS_2020-03-23_17_30_04_3.dat','nchannels',41,'channels',1,'frequency',2000,'duration',17000)
d4 = LoadBinary('AC75a-5_DOB_072519_TS_2020-03-23_17_30_04_4.dat','nchannels',41,'channels',1,'frequency',2000,'duration',17000)

tsData = [d1' ;d2';d3';d4'];

clear d1 d2 d3 d4
tsData = double(tsData);

%%
ix = 1;
A = nan(20000,16);


for i = 1:1000:size(tsData,2)
    At = ltv_adjacency(tsData(:,i:i+4000)');
    A(ix,:) = At(:);
    ix = ix+1
end

xx = tsne(A);

%%

seizFil = 'E:\data\IHKA\Annotations\AC75a-5 DOB 072519_TS_2020-03-23_17_30_04.txt';
TSname1 = 'Seizure starts'; % TS name in file.
TSname2 = 'ends';
TSdata = readtable(seizFil);
TSdata = table2cell(TSdata);
seizure_start = cell2mat(TSdata(cellfun(@any,regexp(TSdata(:,6),TSname1)),4));
seizure_end = cell2mat(TSdata(cellfun(@any,regexp(TSdata(:,6),TSname2)),4));



%%


ts = (1:size(tsData,2))/2000 + 7000;
time2seizure = ts;
kp = true(size(ts));
for i = 1:length(seizure_start)
time2seizure(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
time2seizure(kp& ts>seizure_start(i) & ts<seizure_end(i)) = .5;
time2seizure(kp& ts>seizure_end(i) & ts<seizure_end(i)+600) = 1.5;
time2seizure(kp& ts>seizure_end(i) & ts<seizure_end(i)+600) = 1.5;
kp(ts<seizure_start(i) | ...
(ts>seizure_start(i) & ts<seizure_end(i)) | ...
(ts>seizure_end(i) & ts<seizure_end(i)+600)) = false;
end
time2seizure = time2seizure(1:1000:end);

%%
close all
figure
edg = -20000:100:0;
[~,bb] = histc(time2seizure,edg);
col = linspecer(201,'jet');
for i = 1:size(xx,1)
    if bb(i)>0
    plot(xx(i,1),xx(i,2),'.','color',col{bb(i)},'markersize',10)
    end
hold on
i
end
colormap(cell2mat(col))
colorbar('Ticks',(0:10:100)/100,'TickLabels',edg(1:20:end))

