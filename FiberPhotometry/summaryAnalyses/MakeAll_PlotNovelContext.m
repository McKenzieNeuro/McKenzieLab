load('R:\DANEHippocampalResponse\Workspaces\novelContext.mat')


%%

dayNum = cellfun(@(a) cell2mat(a(:,2)),context1,'uni',0);

clear DA_PETH_NC NE_PETH_NC DA_PETH_HC NE_PETH_HC
for i = 1:11
    
DA_PETH_NC(i,:) = nanmean(cell2mat(cellfun(@(a,b) a(b==i,:),PETH(kp_DA),dayNum(kp_DA),'uni',0)'));
NE_PETH_NC(i,:) = nanmean(cell2mat(cellfun(@(a,b) a(b==i,:),PETH(kp_NE),dayNum(kp_NE),'uni',0)'));

DA_N_NC(i) = sum(cell2mat(cellfun(@(a) sum(a==i),dayNum(kp_DA),'uni',0)'));
NE_N_NC(i) = sum(cell2mat(cellfun(@(a) sum(a==i),dayNum(kp_NE),'uni',0)'));


end

%%
figure
imagesc(ts_PETH,[],DA_PETH)

figure
imagesc(ts_PETH,[],NE_PETH)
%%