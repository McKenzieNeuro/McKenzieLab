clear all
NE = load('R:\DANEHippocampalResponse\Workspaces\novelContext.mat','usubjs','subjID','inRear','beta_NE','tau_hat_NE','beta_HC','tau_hat_HC','beta_NE_days','tau_hat_NE_days');
SOR = load('R:\DANEHippocampalResponse\Workspaces\SOR.mat','usubjs','subjID','beta_SOR','tau_hat');
LT = load('R:\DANEHippocampalResponse\Workspaces\LinearTrack.mat','usubjs','subjID','beta_LT','tau_hat_LT','beta_HC','tau_hat_HC');
%%
all_subj_NE = NE.usubjs(unique(NE.subjID(~(isnan(NE.inRear)))));
all_subj_SOR = SOR.usubjs(unique(SOR.subjID));
all_subj_LT = LT.usubjs(unique(LT.subjID));

%change the name of the subject
all_subj_NE = strrep(all_subj_NE,'NE2m3' ,'NE2h2');


%%
[all_subj,bb,cc] = unique([all_subj_SOR;all_subj_NE;all_subj_NE;all_subj_LT;all_subj_LT]);
beta_all = [SOR.beta_SOR LT.beta_LT LT.beta_HC NE.beta_NE NE.beta_HC]';
tau_hat_all = [SOR.tau_hat LT.tau_hat_LT LT.tau_hat_HC NE.tau_hat_NE NE.tau_hat_HC]';

cond = [ones(length(SOR.beta_SOR),1);2*ones(length( LT.beta_LT),1);3*ones(length( LT.beta_HC),1);4*ones(length( NE.beta_NE),1);5*ones(length( NE.beta_HC),1)];

%%
figure 
hold on
plot(cc(cond==1),tau_hat_all(cond==1),'x')
plot(cc(cond==2),tau_hat_all(cond==2),'o')
plot(cc(cond==3),tau_hat_all(cond==3),'d')
plot(cc(cond==3),tau_hat_all(cond==3),'d')
plot(cc(cond==3),tau_hat_all(cond==3),'d')
set(gca,'xtick',1:8,'xticklabel',all_subj)
legend('SOR','NE','LT')


%%

dec_SOR = beta_all(cond==1).*exp(tau_hat_all(cond==1).*(1:500));
dec_NE = beta_all(cond==2).*exp(tau_hat_all(cond==2).*(1:500));
dec_LT = beta_all(cond==4).*exp(tau_hat_all(cond==4).*(1:500));

%%
figure 
hold on
p(1) = plotMeanSEM(1:500,dec_SOR,'k');
p(2) = plotMeanSEM(1:500,dec_NE,'r');
p(3) = plotMeanSEM(1:500,dec_LT,'b');
legend(p,'SOR','Novel Env','Linear Track')
xlabel('Time (s) from Obj/Con/Track')
ylabel('NE')
%%
close all

for i = [1 2 3 4 5]
    
    errorbar(nanmedian(beta_all(cond==i)),nanmedian(tau_hat_all(cond==i)),SEM((tau_hat_all(cond==i))),SEM((tau_hat_all(cond==i))),SEM(beta_all(cond==i)),SEM(beta_all(cond==i)))
    hold on
end

legend({'SOR','LT-track','LT-HC','NE','home-NE'})

xlabel('\beta')
ylabel('\tau')

%%

figure
col = linspecer(10,'jet');
clear NE_hat
for i = 1:10
    NE_hat(i,:,:) = [NE.beta_NE_days(:,i).*exp(NE.tau_hat_NE_days(:,i).*(1:500))]';
    errorbar(nanmean(NE.beta_NE_days(:,i)),nanmean(NE.tau_hat_NE_days(:,i)),SEM((NE.tau_hat_NE_days(:,i))),SEM((NE.tau_hat_NE_days(:,i))),SEM(NE.beta_NE_days(:,i)),SEM(NE.beta_NE_days(:,i)),'color',col{i})
    hold on
end
ylim([-.12 0])


%%
figure
col = linspecer(10,'jet');
for i = 1:10
    plot(nanmean(NE_hat(i,:,:),3),'color',col{i})
    hold on
end