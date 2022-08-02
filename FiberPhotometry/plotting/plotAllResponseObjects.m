%% change for how the R drive is mounted on your computer
masterDir = 'R:\';
masterDir =  'R:\McKenzieLab\DANEHippocampalResponse';

%% Objects = toys
fils = getAllExtFiles(masterDir,'mat',1);


kp_LT = (cellfun(@any,regexp(fils,'noveltySignal'))) &(cellfun(@any,regexpi(fils,'linear'))); 
kp_SOR = (cellfun(@any,regexp(fils,'novelObject'))) &(cellfun(@any,regexpi(fils,'SOR'))); 


fils = fils(kp_LT);
[dirs,b]  =fileparts(fils);
% 
%       7/27/2022
% SOR
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h1\SOR\DA2h1-211004-112717'                    }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h3\SOR\DA2h3-211007-120617'                    }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\SOR\DA2h4\DA2h4-210921-154048'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DASCD1\SOR\DACSD1_220527\TDT\Subject1-220527-160113'}
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\SOR\NE2h2\NE2h2-211007-125834'}
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\SOR\NE2h4-211213-103915'                    }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\SOR\NE2h6\NE2h6-220511-120129'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\SOR\NE2h7-220322-132227'                    }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\SOR\NE2h8\NE2h8-220321-112022'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\SOR\NE2h9-220419-132755'                    }


% LINEAR TRACK
%     7/29/2022
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h3\LinearTrack\DA2h3-220623-084843\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h3\LinearTrack\DA2h3-220630-101007\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220623-084843\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\LinearTrack\DA2h4-220630-104634\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220623-130352\noveltySignal.mat'}
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220630-112150\noveltySignal.mat'}
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220622-092231\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220623-112508\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220630-125426\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220720-090530\noveltySignal.mat'             }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220721-074443\noveltySignal.mat'             }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220720-093937\noveltySignal.mat'             }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220721-082011\noveltySignal.mat'             }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220622-100042\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220623-121635\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220630-115938\noveltySignal.mat'              }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220720-102405\noveltySignal.mat'             }
%     {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220721-093903\noveltySignal.mat'             }
%%


 % for SOR
%labels = {'Obj1','Obj2','Obj3','Obj4','Obj5'};
%evFile = 'novelObject.mat';

% % for linear track
 labels = {'novelty_on'}; 
 evFile = 'noveltySignal.mat';

%labels = {'Obj1on','Obj2on','Obj3on','Obj4on','Obj5on'};
for i = 1:length(dirs)
    
    cd(dirs{i})
    data = TDTbin2mat(pwd);
    
    if isfield(data.streams,'x405A')
        streams = {'x405A','x465A'};
        isosbestic = 'x405A';
    elseif isfield(data.streams,'x450D')
        streams = {'x450D','x500D'};
        isosbestic = 'x450D';
        
    else
         streams = {'Dv4D','Dv5D'};
        isosbestic = 'Dv4D';
        
    end
    
    [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,evFile,labels,'plotIntervals',[30 30],'streams',streams,'isosbestic',isosbestic);
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
plot(ts_data,nanconvn(signal_DFoF,kernel'),'k')
hold on
ylim([-1 5])
set(gca,'box','off')
xlabel ('Time (s)')
ylabel ('Z-Score')


for i = 1:length(labels)
   plot([ev_tims{i} ev_tims{i}]',[ones(length(ev_tims{i}),1) 2*ones(length(ev_tims{i}),1)]','color',col{i})
end

%%
%get index at time -1
[a,ix_0] = bestmatch(0,ts_PETH);

%get index at time 1
[a,ix_30] = bestmatch(5,ts_PETH);

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

u_obj = cellfun(@(a) double(nanPad(a,5)),u_obj,'uni',0);

%%
col = linspecer(5,'heat');
close all
for j  =  1:size(u_obj,1)
    figure
for i = 1:1
plot(u_obj{j,i}(1),u_obj{j,i}(2),'o','color',col{i})
shg
xlim([-.25 1])
ylim([-.25 1])
hold on
plot(-.15:.01:1,-.15:.01:1)


%plot(u_obj{j,i}(1),u_obj{j,i}(3),'x','color',col{i})
%plot(u_obj{j,i}(1),u_obj{j,i}(4),'d','color',col{i})


xlabel('GRAB 1st encounter')
ylabel('GRAB Nth encounter')

set(gca,'box','off')
end

end


%%
col = flipud(linspecer(5,'jet'));
Nth_to_compare = 2;
% plot mean for DA/NE


kp_DA = cellfun(@any,regexp(dirs,'DA2')) | cellfun(@any,regexp(dirs,'DACSD1'));
kp_NE = cellfun(@any,regexp(dirs,'NE2'));
for ii = 1:2
    
    if ii==1
        kp = kp_NE;
    else
        kp = kp_DA;
    end
    first = cellfun(@(a) a(1), u_obj(kp,:));
    second =  cellfun(@(a) a(Nth_to_compare), u_obj(kp,:));
    u_obj1_DA = nanmean(first);
    u_obj2_DA = nanmean(second);
    
    
    s_obj1_DA = SEM(cellfun(@(a) a(1), u_obj(kp,:)));
    s_obj2_DA = SEM(cellfun(@(a) a(Nth_to_compare), u_obj(kp,:)));
    close all
    h=figure;
    
    for i = 1:5
        
        hold on
        errorbarxy(u_obj1_DA(i),u_obj2_DA(i),s_obj1_DA(i),s_obj1_DA(i),s_obj2_DA(i),s_obj2_DA(i),{'k.-', col{i}, col{i}})
        plot(u_obj1_DA(i),u_obj2_DA(i),'.','color',col{i})
    end
    shg
    xlim([-.5 1.5])
    ylim([-.5 1.5])
    hold on
    plot(-.5:.01:1.5,-.5:.01:1.5)
    
    if ii==1
           xlabel('NE 1st encounter')
    ylabel('NE 2nd encounter')
    else
        
    xlabel('DA 1st encounter')
    ylabel('DA 2nd encounter')
    end
    set(gca,'fontsize',20)
    if ii==2
        saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\DA_SOR_1_vs_2.eps','epsc')
    saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\DA_SOR_1_vs_2.png','png')
    saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\DA_SOR_1_vs_2.fig','fig')
else
    saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\NE_SOR_1_vs_2.eps','epsc')
    saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\NE_SOR_1_vs_2.png','png')
    saveas(h,'R:\McKenzieLab\DANEHippocampalResponse\Figures\NE_SOR_1_vs_2.fig','fig')
end
end
%%
N_obj  = size(u_obj,2);
data = [first(:) ; second(:)];
g1 = [ones(length(first(:)),1); 2*ones(length(first(:)),1)];
g2 = [repmat([1:size(first,1)]',size(first,2),1);repmat([1:size(first,1)]',size(first,2),1)];
g3 = [upSample(1:N_obj,size(first,1));upSample(1:N_obj,size(first,1))];
group = {g1,g2,g3};
[p,tbl,stats,terms] = anovan(data,group,'model','interaction','random',2);



%%
 kp = kp_DA;
    first = cellfun(@(a) a(1), u_obj(kp,:));
    second =  cellfun(@(a) a(Nth_to_compare), u_obj(kp,:));
% linear track
N_obj  = size(u_obj,2);
data = [first(:) ; second(:)];
g1 = [ones(length(first(:)),1); 2*ones(length(first(:)),1)];
g2 = [repmat([1:size(first,1)]',size(first,2),1);repmat([1:size(first,1)]',size(first,2),1)];

group = {g1,g2};
[p,tbl,stats,terms] = anovan(data,group,'model','interaction','random',2);


%%
close all
figure
ax = tight_subplot(5,1);
for i = 1:5
axes(ax(i))
imagesc(ts_PETH,[],newObj1{i},[-1.5 1])

xlabel('Time from Exposure(s)')
set(gca,'box','off')

c = colorbar;
w = c.LineWidth;
c.LineWidth = 0.5;
end

%%
nSamples = size(newObj1{1,1},2);
newObj2 = cellfun(@(a) nanPad(nanPad(a',5)',nSamples),newObj1,'uni',0);

    tmp = nan(5,nSamples,5,9);
%loop over subject
for i = 1:size(newObj2,1)

    %loop over object
    for j = 1:size(newObj2,2)
        
        %loop over sample
        
        for k = 1:5
    tmp(j,:,k,i)  = newObj2{i,j}(k,:);
    
        end
    end
    
end

% tmp = objs x time x sample # x subj

noObj = squeeze(nanmean(tmp,1));

%%
close all
figure

ax = tight_subplot(1,5);


for i = 1:5
    
axes(ax(i))
hold on
plotMeanSEM(ts_PETH,squeeze(noObj(:,i,kp_NE))','k')



plot([0 0],[-.5 2.5],'k')
plot([ts_PETH(1) ts_PETH(end)],[0 0 ],'k')
ylim([-.5 1.5])
xlim([-30 30])
xlabel('Time from object sample (s)')
set(gca,'xtick',[-30:15:30],'xticklabel',-30:15:30,'fontsize',20)
if i ==1
ylabel('NE response')
end

switch i
    case 1
        exp = '1st';
        
    case 2
        exp = '2nd';
        
    case 3
        exp = '3rd';
        
    case 4
        exp = '4th';
        
    case 5
        exp = '5th';
end

title([exp ' sample'])
end