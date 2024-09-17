load('allfeat.mat')
%%
for i = 1:6
    
    dat{i} = cell2mat(arrayfun(@(a) a.dat{i},v,'uni',0)');
    sesID{i} = cell2mat(arrayfun(@(a) a.sesID{i},v,'uni',0)');
end
%%
d = cell2mat(dat');
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat(1:6),num2cell(1:length(dat(1:6))),'uni',0)');

ses = cell2mat(sesID');

%%

ok = tsne(d);

%%
close all
figure
k = gaussian2Dfilter([ 100 100],[10 10]);
for i = 1:6
    subplot(2,3,i)
    bin = histcn(ok(group==i,:),-100:100,-100:100);
    bin = nanconvn(bin,k);
   imagesc(bin) 
end