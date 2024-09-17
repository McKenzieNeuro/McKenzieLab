clear all
NE = load('R:\DANEHippocampalResponse\Workspaces\novelContext_doubleExp.mat','usubjs','subjID','inRear','b_all_NE','b_all_HC','b_all_NE_day');
SOR = load('R:\DANEHippocampalResponse\Workspaces\SOR_doubleExp.mat','usubjs','subjID','b_all_SOR');
LT = load('R:\DANEHippocampalResponse\Workspaces\LinearTrack_doubleExp.mat','usubjs','subjID','b_all_LT','b_all_hc');
%%
all_subj_NE = NE.usubjs(unique(NE.subjID));
all_subj_SOR = SOR.usubjs(unique(SOR.subjID));
all_subj_LT = LT.usubjs(unique(LT.subjID));

%change the name of the subject
%all_subj_NE = strrep(all_subj_NE,'NE2m3' ,'NE2h2');
%all_subj_NE = unique(all_subj_NE);

%%
[all_subj,bb,cc] = unique([all_subj_SOR;all_subj_LT;all_subj_LT;all_subj_NE;all_subj_NE]);

beta_all = [SOR.b_all_SOR(:,1); LT.b_all_LT(:,1); LT.b_all_hc(:,1); NE.b_all_NE(:,1) ; NE.b_all_HC(:,1);NE.b_all_NE(:,9)];
tau_hat_all = [SOR.b_all_SOR(:,2); LT.b_all_LT(:,2); LT.b_all_hc(:,2); NE.b_all_NE(:,2);  NE.b_all_HC(:,2);NE.b_all_NE(:,10)];

beta_all_neg = [SOR.b_all_SOR(:,3); LT.b_all_LT(:,3); LT.b_all_hc(:,3); NE.b_all_NE(:,3) ; NE.b_all_HC(:,3);NE.b_all_NE(:,10)];
tau_hat_all_neg = [SOR.b_all_SOR(:,4); LT.b_all_LT(:,4); LT.b_all_hc(:,4); NE.b_all_NE(:,4);  NE.b_all_HC(:,4);NE.b_all_NE(:,12)];


cond = [ones(size(SOR.b_all_SOR,1),1);2*ones(size( LT.b_all_LT,1),1);3*ones(size( LT.b_all_hc,1),1);4*ones(size( NE.b_all_NE,1),1);5*ones(size( NE.b_all_HC,1),1);6*ones(size( NE.b_all_NE,1),1)];

%%
figure 
hold on
plot(cc(cond==1),tau_hat_all(cond==1),'x')
plot(cc(cond==2),tau_hat_all(cond==2),'o')
plot(cc(cond==3),tau_hat_all(cond==3),'d')
plot(cc(cond==4),tau_hat_all(cond==4),'.')
plot(cc(cond==5),tau_hat_all(cond==5),'p')
set(gca,'xtick',1:8,'xticklabel',all_subj)
legend('SOR','NE','LT')


%%

dec_SOR = beta_all(cond==1).*exp(tau_hat_all(cond==1).*(1:500)) + beta_all_neg(cond==1).*exp(tau_hat_all_neg(cond==1).*(1:500));
dec_NE = beta_all(cond==2).*exp(tau_hat_all(cond==2).*(1:500)) + beta_all_neg(cond==2).*exp(tau_hat_all_neg(cond==2).*(1:500));
dec_LT = beta_all(cond==4).*exp(tau_hat_all(cond==4).*(1:500)) + beta_all_neg(cond==4).*exp(tau_hat_all_neg(cond==4).*(1:500));

dec_HC_NE = beta_all(cond==3).*exp(tau_hat_all(cond==3).*(1:500))  + beta_all_neg(cond==3).*exp(tau_hat_all_neg(cond==3).*(1:500));
dec_HC_LT = beta_all(cond==5).*exp(tau_hat_all(cond==5).*(1:500))  + beta_all_neg(cond==5).*exp(tau_hat_all_neg(cond==5).*(1:500));


dec_rear = beta_all(cond==6).*exp(tau_hat_all(cond==6).*(1:500)) + beta_all_neg(cond==6).*exp(tau_hat_all_neg(cond==6).*(1:500));

%%

beta_all_neg1 =accumarray([cc cond],beta_all_neg,[],@sum,nan);

beta_all1 =accumarray([cc cond],beta_all,[],@sum,nan);


tau_hat_all1 =accumarray([cc cond],tau_hat_all,[],@sum,nan);

tau_hat_all_neg1 =accumarray([cc cond],tau_hat_all_neg,[],@sum,nan);


%%
figure 
hold on
p(1) = plotMeanSEM(1:500,dec_SOR,'k');
p(2) = plotMeanSEM(1:500,dec_NE,'r');
p(3) = plotMeanSEM(1:500,dec_LT,'b');
p(4) = plotMeanSEM(1:500,dec_rear,'y');

legend(p,'SOR','Novel Env','Linear Track','rear')
xlabel('Time (s) from Obj/Con/Track')
ylabel('NE')
%%
figure 
hold on
p(1) = plotMeanSEM(1:500,dec_HC_NE,'r');
p(2) = plotMeanSEM(1:500,dec_HC_LT,'k');
%%
close all

for i = [1 2 3 4 5 6]
    
    errorbar(nanmedian(beta_all(cond==i)),nanmedian(tau_hat_all(cond==i)),SEM((tau_hat_all(cond==i))),SEM((tau_hat_all(cond==i))),SEM(beta_all(cond==i)),SEM(beta_all(cond==i)))
    hold on
end

legend({'SOR','LT-track','LT-HC','NE','home-NE','rear'})

xlabel('\beta')
ylabel('\tau')

%%
close all
figure
col = linspecer(7,'gray');
for i = [3 5]
  
    errorbar(nanmedian(beta_all_neg(cond==i)),nanmedian(tau_hat_all_neg(cond==i)),SEM((tau_hat_all_neg(cond==i))),SEM((tau_hat_all_neg(cond==i))),SEM(beta_all_neg(cond==i)),SEM(beta_all_neg(cond==i)),'color',col{i})
    hold on
end

legend({'SOR','LT-track','LT-HC','NE','home-NE'})

xlabel('\beta')
ylabel('\tau')


%%

figure
col = linspecer(10,'jet');
clear NE_hat
for i = 1:8
    for j = 1:10
    NE_hat(i,:,j) = NE.b_all_NE_day{j}(i,1).*exp(NE.b_all_NE_day{j}(i,2).*(1:500)) + ...
        NE.b_all_NE_day{j}(i,3).*exp(NE.b_all_NE_day{j}(i,4).*(1:500));
    %hold on
    end
end
%%
close all
i=3;
    errorbar(nanmedian(beta_all(cond==i)),nanmedian(tau_hat_all(cond==i)),SEM((tau_hat_all(cond==i))),SEM((tau_hat_all(cond==i))),SEM(beta_all(cond==i)),SEM(beta_all(cond==i)),'color','k')
hold on
for j = 1:10
     errorbar(nanmedian(NE.b_all_NE_day{j}(:,1)),nanmedian(NE.b_all_NE_day{j}(:,2)),SEM(NE.b_all_NE_day{j}(:,2)), ...
         SEM(NE.b_all_NE_day{j}(:,2)),SEM(NE.b_all_NE_day{j}(:,1)), ...
         SEM(NE.b_all_NE_day{j}(:,1)),'color',col{j})
   
     hold on
%ylim([-.06 0])
end
lab = cellfun(@(a) ['day ' num2str(a)],num2cell(1:10),'UniformOutput',false);
lab = [{'Home cage'};lab'];
legend(lab)
xlabel('\beta')
ylabel('\tau')

%%

col = linspecer(10,'jet');

close all
i=3;
    errorbar(nanmedian(beta_all_neg(cond==i)),nanmedian(tau_hat_all_neg(cond==i)),...
        SEM((tau_hat_all_neg(cond==i))),SEM((tau_hat_all_neg(cond==i))),...
        SEM(beta_all_neg(cond==i)),SEM(beta_all_neg(cond==i)),'color','k')
hold on
for j = 1:10
     errorbar(nanmedian(NE.b_all_NE_day{j}(:,3)),nanmedian(NE.b_all_NE_day{j}(:,4)),SEM(NE.b_all_NE_day{j}(:,4)), ...
         SEM(NE.b_all_NE_day{j}(:,4)),SEM(NE.b_all_NE_day{j}(:,3)), ...
         SEM(NE.b_all_NE_day{j}(:,3)),'color',col{j})
   
     hold on
%ylim([-.06 0])
end
lab = cellfun(@(a) ['day ' num2str(a)],num2cell(1:10),'UniformOutput',false);
lab = [{'Home cage'};lab'];
legend(lab)
xlabel('\beta')
ylabel('\tau')

%%
hold on
col = linspecer(10,'jet');

for j = 1:10
     errorbar(nanmean(NE.b_all_NE_day{j}(:,3)),nanmean(NE.b_all_NE_day{j}(:,4)), ...
         SEM(NE.b_all_NE_day{j}(:,4)), ...
         SEM(NE.b_all_NE_day{j}(:,4)), ...
         SEM(NE.b_all_NE_day{j}(:,3)), ...
         SEM(NE.b_all_NE_day{j}(:,3)),'color',col{j})
   
     hold on
ylim([-.12 0])
end


%%
figure
col = linspecer(10,'jet');
for i = 1:10
    plot(nanmean(NE_hat(:,:,i),1),'color',col{i})
    hold on
end
figure
for i = 1:10
plot(0:1:500,avghist(timeFrmEntry(kp_Novel==1 &   dayNum==i),neural(kp_Novel==1 &   dayNum==i),0:1:500),'color',col{i})
hold on
end