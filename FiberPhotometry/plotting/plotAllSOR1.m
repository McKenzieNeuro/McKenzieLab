% plot DA SOR
clear all
close all
dirs = [ ...

{'R:\McKenzieLab\DANEHippocampalResponse\DA2h3\SOR\DA2h3-211007-120617'} ;...
{'R:\McKenzieLab\DANEHippocampalResponse\DA2h4\SOR\DA2h4\DA2h4-210921-154048'} ;...
{'R:\McKenzieLab\ASommer\GNemer\Data\Novel Object\Single Day\DA2h9\DA2h9-230215-110140'}; ...
{'R:\McKenzieLab\DANEHippocampalResponse\DASCD1\SOR\DACSD1_220527\TDT\Subject1-220527-160113'} ; ...
{'R:\McKenzieLab\ASommer\GNemer\Data\Novel Object\Multi Day\DA3h4'}]; % has many sub directories


% dirs = [ ...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\SOR\NE2h4-211213-103915'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h4\SOR\NE2h4-211213-103915'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h6\SOR\NE2h6-220511-120129'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h7\SOR\NE2h7-220322-132227'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h8\SOR\NE2h8\NE2h8-220321-112022'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h9\SOR\NE2h9-220419-132755'};...
%    {'R:\McKenzieLab\DANEHippocampalResponse\NE2h10\SOR\NE2h10-230215-115346'}];

%%

clear newObj

%labels = {'Obj1on','Obj2on','Obj3on','Obj4on','Obj5on'};
for i = 1:length(dirs)
    
    if any(cellfun(@any,getAllExtFiles(dirs{i},'tev',0)))
        goodDir{1} = dirs{i};
        
    else
        fils = getAllExtFiles( dirs{i},'tev',1);
        goodDir = fileparts(fils);
    end
    
    
    
    for k = 1:length(goodDir)
        clear evFil
        cd(goodDir{k})
        data = TDTbin2mat(goodDir{k});
        
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
        
        if exist('novelObjectIntro.mat')
            evFil = 'novelObjectIntro.mat';
        elseif exist('ContextTransitionRevised.mat')
            evFil = 'ContextTransitionRevised.mat';
        elseif exist('objects1.mat')
            evFil = 'objects1.mat';
        end
        


        v = load(evFil);
        evType = sort(strrep(unique(v.data(:,1)),' ',''));
        evType = evType(~(cellfun(@any,regexp(evType,'off'))) & ~(cellfun(@any,regexp(evType,'empty'))));
        
        
        [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(goodDir{k},evFil,evType,'plotIntervals',[30 120],'streams',streams,'isosbestic',isosbestic);
        for j = 1:length(evType)
            newObj{i,j,k} = signal_DFoF(ix{j});
        end
        
        
        
    end
end

%%

k  = gaussian2Dfilter([1 10000 ],1250);

newObj1 = cellfun(@(a) nanconvn(a,k),newObj,'uni',0);
clear o
for  i =1:5
    
    o(:,:,i) = cell2mat(newObj1(:,i,1));
    
end

%%
col = flipud(linspecer(5,'lines'));
set(0, 'DefaultFigureRenderer', 'painters');
outdir = 'R:\ASommer\GNemer\figures\SOR';
figure
ax = tight_subplot(1,5);
for i = 1:5
    axes(ax(i))
        hold on
    plot([-50 250],[0 0],'k')
    plot([0 0],[-.5 1.5],'k')
    plotMeanSEM(ts_PETH(1:10:end),o(:,1:10:end,i),col{i})
    ylim([-.4 1.4])
    xlim([-50 120])
    outfil = [outdir filesep 'obj' num2str(i) ];
    % titlestr = sprintf('Obj %d',i);
    % title(titlestr,'FontSize',14)
    % saveas(gcf,outfil,'epsc')
    % close all

end



%%
% 
% o = squeeze(cell2mat(newObj1(5,1,:)))';
% for i = 1:5
%     figure
% plot(ts_PETH(1:10:end),o(i,1:10:end))
%   ylim([-.4 2.4])
%   outfil = [outdir filesep 'day' num2str(i) ];
%     % saveas(gcf,outfil,'epsc')
%     close all
% end
