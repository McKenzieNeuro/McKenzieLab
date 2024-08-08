
% directories with no novelty, no double

directoryNovel = [ ...
    {' R:\McKenzieLab\DANEHippocampalResponse\DA2h3\LinearTrack\DA2h3-220623-084843'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\DA2h3\LinearTrack\DA2h3-220630-101007'} ; ...
    
    {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220623-084843'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220630-104634'} ; ...
    
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220623-130352'} ; ... %(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220630-112150'} ; ...%(done)
    
   % {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220622-092231'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220623-112508'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220630-125426'} ; ...%(done)
    
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220622-100042'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220623-121635'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220630-115938'} ; ...%(done)
    
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220720-090530'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220721-074443'} ; ...%(done)
    
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220720-093937'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220721-082011'} ; ....%(done)
    
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220720-102405'} ; ...%(done)
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220721-093903'} ; ...%(done)
    ];



directoryControl = [ ...
    {'  R:\McKenzieLab\DANEHippocampalResponse\DA2h3\LinearTrack\DA2h3-220701-110500'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220624-130830'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220701-113832'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220624-115525'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220701-090608'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220624-102557'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220701-094335'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220624-110628'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220701-101956'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220722-091944'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220722-095135'} ; ...
    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220722-102506'} ; ...
    
    ];



%nSessions = length(directory);
%%
%clear left right ts_PETH
for i = 1:length(directoryNovel)
    if exist(directoryNovel{i})
   [home_resp1{i},home_resp2{i},track_respNO{i},track_respLT1{i},ts_PETH] = plotTrackResponseNovel(directoryNovel{i});
    i
    end
end
%save('double_reward_linearTrack.mat','-v7.3')

%%

for i = 1:length(directoryControl)
    if exist(directoryControl{i})
   [home_resp1_C{i},home_resp2_C{i},track_respNO_C{i},track_respLT1_C{i},ts_PETH] = plotTrackResponseNovel(directoryControl{i});
    i
    end
end

%%

kp_DA = cellfun(@any,regexp(directoryNovel,'DA2'));
kp_NE = cellfun(@any,regexp(directoryNovel,'NE2'));
kp =~cellfun(@isempty,home_resp1) & kp_DA';
hc1 = cell2mat(home_resp1(kp)');
hc2 = cell2mat(home_resp2(kp)');
NO = cell2mat(track_respNO(kp)');
LT1 = cell2mat(track_respLT1(kp)');

close all
figure
hold on


plotMeanSEM(ts_PETH{4}(1:100:end),LT1(:,1:100:end),'k')
plotMeanSEM(ts_PETH{1}(1:100:end)+275,hc1(:,1:100:end),'k')
plotMeanSEM(ts_PETH{2}(1:100:end)+1100,hc2(:,1:100:end),'k')
plotMeanSEM(ts_PETH{3}(1:100:end)+400,NO(:,1:100:end),'k')

ts1 = [ts_PETH{4}(1):100:ts_PETH{4}(end)];
ts2 = [-30 0 30]+200;
ts3 = [ 0 100 200 300]+275;
ts4 = [ts_PETH{2}(1):100:ts_PETH{2}(end)]+1100;
ts = round([ts1 ts2 ts3 ts4]);

ts1l = [ts_PETH{4}(1):100:ts_PETH{4}(end)];
ts2l = [-30 0 30 ];
ts3l = [ 0 100 200 300];
ts4l = [ts_PETH{2}(1):100:ts_PETH{2}(end)];
tsl = round([ts1l ts2l ts3l ts4l]);


set(gca,'xtick',ts,'xticklabel',tsl)




kp_DA = cellfun(@any,regexp(directoryControl,'DA2'));
kp_NE = cellfun(@any,regexp(directoryControl,'NE2'));
kp =~cellfun(@isempty,home_resp1_C) & kp_DA';
hc1 = cell2mat(home_resp1_C(kp)');
hc2 = cell2mat(home_resp2_C(kp)');
NO = cell2mat(track_respNO_C(kp)');
LT1 = cell2mat(track_respLT1_C(kp)');

hold on


plotMeanSEM(ts_PETH{4}(1:100:end),LT1(:,1:100:end),'r')
plotMeanSEM(ts_PETH{1}(1:100:end)+275,hc1(:,1:100:end),'r')
plotMeanSEM(ts_PETH{2}(1:100:end)+1100,hc2(:,1:100:end),'r')
plotMeanSEM(ts_PETH{3}(1:100:end)+400,NO(:,1:100:end),'r')

ts1 = [ts_PETH{4}(1):100:ts_PETH{4}(end)];
ts2 = [-30 0 30]+275;
ts3 = [ 0 100 200 300]+400;
ts4 = [ts_PETH{2}(1):100:ts_PETH{2}(end)]+1100;
ts = round([ts1 ts2 ts3 ts4]);

ts1l = [ts_PETH{4}(1):100:ts_PETH{4}(end)];
ts2l = [-30 0 30 ];
ts3l = [ 0 100 200 300];
ts4l = [ts_PETH{2}(1):100:ts_PETH{2}(end)];
tsl = round([ts1l ts2l ts3l ts4l]);


set(gca,'xtick',ts,'xticklabel',tsl)

%%

kp_DA = cellfun(@any,regexp(directory,'DA2'));
kp_NE = cellfun(@any,regexp(directory,'NE2'));% & cellfun(@any,regexp(directory,'NE2m3'));  





for i = 2
    
    if  i ==1
        kp = kp_DA;
    else
        kp = kp_NE;
        
    end
    singles1 = cellfun(@double,singles(kp),'uni',0);
    doubles1 = cellfun(@double,doubles(kp),'uni',0);
    other_side1 = cellfun(@double,other_side(kp),'uni',0);
    
    singles1 = cell2mat(singles1');
    doubles1 = cell2mat(doubles1');
    other_side1 = cell2mat(other_side1');
    h = figure;
    ax  = tight_subplot(3,1);
    for j = 1:3
        axes(ax(j))
        switch j
            case 1
                data = singles1;
            case 2
                data = doubles1 ;
            case 3
                data = other_side1;
        end
        
        
        %  data  =[left1;right1];
          data = nanzscore(data,[],2);
        
        plotMeanSEM(ts_PETH,data,'k')
        
   % ylim([-1 2])    
    end
    
    if i ==1
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_DAZ.fig','fig')
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_DAZ.png','png')
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_DAZ.eps','epsc')
        
        
    else
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_NEZ.fig','fig')
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_NEZ.png','png')
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\double_reward_mean_NEZ.eps','epsc')
    end
    
    close all
end

