function alpha=sm_power_repeated_anova(mu,sigma,nsub,pct,ratio_within_across)
%mu  = mean per condition 
nrep = length(pct);
for n = 1:1000
    dat = normrnd( mu , sigma,[nsub 1] ) ;
for i = 2:nrep
    dat = [dat dat(:,1) + normrnd( pct(i) * dat(:,1),(1+pct(i))*sigma/ratio_within_across)];
end
group = repmat(1:nrep,nsub,1);
subj = repmat([1:nsub]',1,nrep);
 [pt,t,stats,terms] = anovan(dat(:),[{group(:)} {subj(:)}],'random',[2],'display',false);
p(n) = pt(1);
end

alpha = nanmean(p<.05);
end