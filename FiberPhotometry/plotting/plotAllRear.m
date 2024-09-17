fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);
fils = fils(contains(fils,'rear.mat'));
%%
k  = gaussian2Dfilter([10000 1],1250);
for i  = 1:length(fils)
    cd(fileparts(fils{i}))
    [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,'rear.mat',{'startRear'});
    signal_DFoF = nanconvn(signal_DFoF,k');
    rear_PETH{i} = signal_DFoF(ix{1});
    i
    
end

%%
idx = regexp(fils,'NE');
mice = cellfun(@(a,b) a(b(2):b(2)+4),fils,idx,'UniformOutput',false);
ok = cell2mat(cellfun(@(a) nanmean(zscore(double(a),[],2)),rear_PETH,'uni',0)');

kp1 = contains(mice,'NE2h4') | contains(mice,'NE2h8') | contains(mice,'NE2h4') |  contains(mice,'NE2h9');
kp2 = ~kp1;
%%
[a,b,c] = unique(mice);
nM=  length(a);
col = linspecer(nM,'jet');
figure
ax  = tight_subplot(2,4);
for i = 1:nM
    axes(ax(i))
    m2 = cell2mat(cellfun(@(a) zscore(double(a),[],2),rear_PETH(c==i),'uni',0)');
    %  m2 = cell2mat(cellfun(@(a) double(a),rear_PETH(c==i),'uni',0)');
   imagesc(ts_PETH,[],m2,[-1 1])
   set(gca,'ydir','normal')
   hold on
   plot(ts_PETH,nanmean(m2)*size(m2,1)/3+size(m2,1)/2,'w','linewidth',4)
title(a{i})
end