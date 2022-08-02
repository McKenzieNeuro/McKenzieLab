function autoplotObj()

topDir = uigetdir();
files = dir(topDir);
dirFlags = [files.isdir];
subdirs = files(dirFlags);
dirs = {subdirs(3:end).name};

%%
clear newObj

labels = {'Obj1','Obj2','Obj3','Obj4','Obj5'};
%labels = {'Obj1on','Obj2on','Obj3on','Obj4on','Obj5on'};
for i = 1:length(dirs)
    
    cd(dirs{i})
    [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,'novelObject.mat',labels,'photoBleachCorrection','exp2','plotIntervals',[10 10],'returnedDataType','corrected');
    for j = 1:length(labels)
        newObj{i,j} = signal_DFoF(ix{j});
    end
    close all 
end

clear ix 

%%
figure
kernel  = gaussian2Dfilter([1000 1],[500 1]);
col = linspecer(length(labels),'jet');
col{1} = [1 0 0];
col{2} = [0 0 1];
col{2} = [0 1 0];
%close all
plot(ts_data,nanconvn(signal_DFoF,kernel'))
hold on
ylim([-1 5])

for i = 1:length(labels)
   plot([ev_tims{i} ev_tims{i}]',[ones(length(ev_tims{i}),1) 2*ones(length(ev_tims{i}),1)]','color',col{i})
end

%%
%get index at time 0
[a,ix_0] = bestmatch(0,ts_PETH);

%get index at time 30
[a,ix_30] = bestmatch(30,ts_PETH);

clear u_obj


newObj1 = newObj;
for i = 1:size(newObj,1)
    
    for j = 1:size(newObj,2)
        
        for k = 1:size(newObj1{i,j},1)
        newObj1{i,j}(k,:) = nanconvn( newObj1{i,j}(k,:) ,kernel');
        
        
        u_obj{i,j}(k) = nanmean( newObj1{i,j}(k,ix_0:ix_30),2);




        end
    end
end

%%
col = linspecer(5,'heat');
close all
for j  =  1:size(u_obj,1)
    figure
for i = 1:5
plot(u_obj{j,i}(1),u_obj{j,i}(2),'o','color',col{i})
shg
xlim([-.25 1])
ylim([-.25 1])
hold on
plot(-.15:.01:1,-.15:.01:1)


plot(u_obj{j,i}(1),u_obj{j,i}(3),'x','color',col{i})
plot(u_obj{j,i}(1),u_obj{j,i}(4),'d','color',col{i})


xlabel('GRAB 1st encounter')
ylabel('GRAB Nth encounter')
end

end
%%


close all
figure
ax = tight_subplot(5,1);
for i = 1:5
axes(ax(i))
imagesc(ts_PETH,[],newObj1{i},[-1.5 1])
title('New Object1')
xlabel('Time from exposure')
end

%%
%close all
figure
for i = 1:5


plot(ts_PETH,nanmean(newObj1{i},1))

hold on



end