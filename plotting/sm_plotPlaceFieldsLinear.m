%get positio
%clear all

warning off
spikes = bz_GetSpikes;


basepath  =pwd;
basename = bz_BasenameFromBasepath(basepath);



load('position_info.mat')



%%
kp = cellfun(@range,pos_inf.len_ep)>50;

pos_inf.in_eps = pos_inf.in_eps(kp);

pos_inf.out_eps = pos_inf.out_eps(kp);

pos_inf.len_ep = pos_inf.len_ep(kp);

pos_inf.ts_ep = pos_inf.ts_ep(kp);


ineps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.in_eps),'uni',0));
outeps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.out_eps),'uni',0));






t= pos_inf.ts;
%kp = InIntervals(sti{1}(:,1),[min(t) max(t)]);
len_in_st =[];
len_in_fin =[];
len_out_st =[];
len_out_fin =[];

interval_out =[];
len = pos_inf.lin_pos;



t= pos_inf.ts;

kp = cellfun(@range,pos_inf.ts_ep)<30;

pos_inf.len_ep = pos_inf.len_ep(kp);
pos_inf.ts_ep = pos_inf.ts_ep(kp);
pos_inf.in_eps = pos_inf.in_eps(kp);
pos_inf.out_eps = pos_inf.out_eps(kp);
ineps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.in_eps),'uni',0));
outeps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.out_eps),'uni',0));





%%




%plot data for each neu1ron
bins = min(pos_inf.lin_pos):range(pos_inf.lin_pos)/100:max(pos_inf.lin_pos);
for i = 	1:length(spikes.times)
    h = figure('visible','on');
    ax = tight_subplot(1,2);
    
    
    axes(ax(1))
    spks1 =  spikes.times{i};
    if length(spks1)<100900
        [status,interval] = InIntervals(spks1,ineps);
        
        n1 = histoc(interval(interval>0),1:size(ineps,1));
        
        
        %     len_ep = mat2cell(len(status),n1);
        
        
        spks = spks1(status);
        [~,b] = histc(spks,t);
        unit(i).linpos = nan(size(spks,1),1);
        unit(i).linpos(b>0) = len(b(b>0));
        
        
        
        unit(i).pos = mat2cell( unit(i).linpos,n1);
        % set(gca,'ytick',1:25,'yticklabel',1:25)
        % pos = cellfun(@(a) inpaint_nans(a),pos,'uni',0);
        
        
        hold on
        %  imagesc(bins,[],cell2mat(cellfun(@(a) histoc(a,bins)',pos_inf.len_ep(pos_inf.in_eps),'uni',0))>0)
        colormap([1 1 1 ; .9 .9 .9])
        
        PipeRaster(unit(i).pos)
        
        
        axis on
        title('inbound: stim')
        xlim([min(bins) max(bins)])
        axes(ax(2))
        if any(outeps)
            [status,interval] = InIntervals(spks1,outeps);
        else
            status = false(length(spks1),1);
        end
        n1 = histoc(interval(interval>0),1:size(outeps,1));
        
        
        
        spks = spks1(status);
        [~,b] = histc(spks,t);
        unit(i).linpos = nan(size(spks,1),1);
        unit(i).linpos(b>0) = len(b(b>0));
        
        
        
        unit(i).pos = mat2cell( unit(i).linpos,n1);
        
        %  pos = cellfun(@(a) inpaint_nans(a),pos,'uni',0);
        
        hold on
        %  imagesc(bins,[],cell2mat(cellfun(@(a) histoc(a,bins)',pos_inf.len_ep(pos_inf.out_eps),'uni',0))>0)
        colormap([1 1 1 ; .9 .9 .9])
        %
         PipeRaster(unit(i).pos)
        xlim([min(bins) max(bins)])
        title('outbound: non-stim')
        
        mtit(['unit ' num2str(i)])
        
        set(h,'Units','Inches');
        
        pos = get(h,'Position');
        
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        
        set(groot,'DefaultFigureRenderer','painters')
        
        % print(h, '-dpsc2', [basename '.ps'], '-append');
        
        uiwait
        close all
    end
end
