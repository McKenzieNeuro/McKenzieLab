v = load('R:\DANEHippocampalResponse\Workspaces\linearTrack_NO_con.mat');
vv = load('R:\DANEHippocampalResponse\Workspaces\linearTrack_NO.mat');
%%
figure

 bar([ nanmean(v.b_all_LT([1:5],6)) nanmean(vv.b_all_LT([1:5],6))],'facecolor','w')
 
 hold on
  plot([ones(5,1) 2*ones(5,1)]',[v.b_all_LT([1:5],6) vv.b_all_LT([1:5],6)]')
 set(gca,'xticklabel',{'Con','Obj'}) 
  %%
  
  close all
plot([v.b_all_LT([1:5],5) vv.b_all_LT([1:5],5)]',[v.b_all_LT([1:5],6) vv.b_all_LT([1:5],6)]')

hold on
plot([v.b_all_LT([1:5],5)]',[ v.b_all_LT([1:5],6)]','o')
