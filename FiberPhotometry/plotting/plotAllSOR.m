
vidScore = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
SOR_NE = vidScore(cellfun(@any,regexpi(vidScore,'object')) & cellfun(@any,regexp(vidScore,'NE2')) & cellfun(@any,regexp(vidScore,'SOR')));

SOR_NE = vidScore(cellfun(@any,regexpi(vidScore,'object')) & (cellfun(@any,regexp(vidScore,'DA2')) | cellfun(@any,regexp(vidScore,'DACSD'))) & cellfun(@any,regexp(vidScore,'SOR')))
%SOR_NE = vidScore(cellfun(@any,regexpi(vidScore,'foodReward')) & cellfun(@any,regexp(vidScore,'DA2')) );



for i = 1:length(SOR_NE)
    
    [TDTDirectory,b]  =fileparts(SOR_NE{i});
    
    v = load(SOR_NE{i});
    evType = sort(strrep(unique(v.data(:,1)),' ',''));
    evType = evType(~(cellfun(@any,regexp(evType,'off'))));
    cd(TDTDirectory)
    
    [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(TDTDirectory,[b '.mat'],evType);
    
    ker  = gaussian2Dfilter([1 10000 ],500);
    signal_DFoF = nanconvn(signal_DFoF,ker);
    
    for j = 1:length(ix)
        
        
        for k = 1:5
            if j==1 && i==1
                obj_resp{k} = nan(length(ix),101726,length(SOR_NE));
                
                
            end
            if size(ix{j},1)>=k
                obj_resp{k}(j,:,i) = signal_DFoF(ix{j}(k,:));
            end
        end
        
        
    end
    i
end


%%
col = flipud(linspecer(5,'jet'));
figure
ax = tight_subplot(1,5);
for i = 1:5
    axes(ax(i))
    hold on
    plot([-50 50],[0 0],'k')
    plot([0 0],[-.5 1.5],'k')
    plotMeanSEM(ts_PETH,squeeze(nanmean(obj_resp{i},1))','k')
    ylim([-.5 1.5])
    xlim([-50 50])
    switch i
        case 1
            title('1st sample')
        case 2
            title('2nd sample')
        case 3
            title('3rd sample')
        case 4
            title('4th sample')
        case 5
            title('5th sample')
    end
    xlabel('time to sample (s)')
    if i==1
        ylabel('NE response')
    end
end