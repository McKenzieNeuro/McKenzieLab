
PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);

%%

clustThresh = 0;
custDist = 10;
combDist = 5;
clustSize = 10;
plotThresh = 0;

clf
n = 2;
searchData = cell(1,n);
for i = 1:n
    
    % Get data from file
    fileName = [PathName,'\Box',num2str(i),'_gridData_TD.mat'];
    gridData = matfile(fileName,'Writable',true);
    
    try
        searchData{1,i} = gridData.search;
    catch
        continue
    end

    numObs = size(gridData.obs,1);
    if numObs > 0
        
        % Get coordinates from parameters
        x = gridData.search(1:numObs,1);
        y = gridData.search(1:numObs,2);
        z = gridData.search(1:numObs,3);
        c = -gridData.obs(:,7);
        cLim = [plotThresh,200];
        
        % Find parameters that have not been tested
        if numObs < size(gridData.search,1)
            notTested = nan(size(gridData.search,1)-numObs,1);
            searchData{1,i}(:,4) = [c;notTested];
        else
            searchData{1,i}(:,4) = c;
        end
        
        % Find entries that are above threshold
        nzIdx = find(c >= clustThresh);
        nzObs = c(nzIdx,1);
        numObs = length(nzObs);
        
        % Calculate distance between all points
        d = pdist([x(nzIdx,:),y(nzIdx,:),z(nzIdx,:)],'euclidean');
        s = squareform(d);
        
        % Group points together within the distance of custDist
        distGroup = zeros(numObs,1);
        for ii = 1:numObs
            pd = find((s(:,ii) > 0)+(s(:,ii) <= custDist) == 2);
            g = [ii,pd'];
            distGroup(ii,1:length(g)) = g;
        end
        
        % Get the mean for each group of points
        mXYZ = zeros(numObs,3);
        for ii = 1:numObs
            idx = nzIdx(distGroup(ii,distGroup(ii,:) > 0),1);
            mXYZ(ii,1) = mean(x(idx,1));
            mXYZ(ii,2) = mean(y(idx,1));
            mXYZ(ii,3) = mean(z(idx,1));
        end
        
        % Calculate distance between all means
        d = pdist([mXYZ(:,1),mXYZ(:,2),mXYZ(:,3)],'euclidean');
        s = squareform(d);
        
        % Group means together within the distance of combDist
        distGroup_clust = zeros(numObs,1);
        for ii = 1:numObs
            pd = find((s(:,ii) <= combDist) == 1);
            g = unique([ii,pd']);
            distGroup_clust(ii,1:length(g)) = g;
        end
     
        % Combine all overlapping groups
        numGroup = 0;
        clust = zeros(numObs,1);
        for ii = 1:numObs
            if clust(ii,1) == 0

                numGroup = numGroup+1;
                clust(ii,1) = numGroup;
                g = distGroup_clust(ii,distGroup_clust(ii,:) > 0);

                for i3 = 1:numObs
                    if i3 ~= ii && clust(i3,1) == 0
                        if any(ismember(distGroup_clust(i3,:),g))
                            g = unique([g,distGroup_clust(i3,distGroup_clust(i3,:) > 0)]);
                            clust(i3,1) = numGroup;
                        end
                    end
                end
            end
        end
        
        % Group all indicies into clusters
        clustIdx = cell(numGroup,1);
        for ii = 1:numObs
            temp = unique(distGroup(distGroup_clust(ii,distGroup_clust(ii,:) > 0),:));
            temp(temp == 0) = [];
            if size(temp,2) > 1
                temp = temp';
            end
            clustIdx{clust(ii,1),1} = [clustIdx{clust(ii,1),1};temp];
        end
        
        % Remove duplicate entries
        for ii = 1:numGroup
            clustIdx{ii,1} = unique(clustIdx{ii,1});
        end
        
        % Shift points that occur in multiple clusters to the largest cluster
        for ii = 1:numGroup
            g = clustIdx{ii,1};
            for i3 = 1:numGroup
                if i3 ~= ii
                    g2 = clustIdx{i3,1};
                    m = ismember(g2,g);
                    for i4 = 1:length(m)
                        if m(i4) == 1
                            [~,im] = min([length(g),length(g2)]);
                            if im == 1
                                clustIdx{ii,1}(i4) = nan;
                            elseif im == 2
                                clustIdx{i3,1}(i4) = nan;
                            end
                        end
                    end
                end
            end
        end
        
        % Find the size of all clusters
        sizeClust = zeros(numGroup,1);
        for ii = 1:numGroup
            clustIdx{ii,1}(isnan(clustIdx{ii,1})) = [];
            clustIdx{ii,1}(clustIdx{ii,1} == 0) = [];
            sizeClust(ii,1) = length(clustIdx{ii,1});
        end
        
        % Remove clusters less than clustSize
        idx = find(sizeClust < clustSize);
        clustIdx(idx) = [];
        sizeClust(idx) = [];

        [sizeClust,idx] = sort(sizeClust,'descend');
        clustIdx = clustIdx(idx,1);

        % Label points with cluster ID
        nzIdx(:,2) = 0;
        for ii = 1:length(clustIdx)
            nzIdx(clustIdx{ii,1},2) = ii;
        end
        
        % Get mean score and parameters for each cluster
        nzIdx(:,3:7) = 0;
        for ii = 1:numGroup

            idx = nzIdx(nzIdx(:,2) == ii,1);
            idx2 = find(nzIdx(:,2) == ii);
            
            score = c(idx,1);
            nzIdx(idx2,3) = round(mean(score),3);
            nzIdx(idx2,4) = (sum(score == 0)/length(score))*100;

            nzIdx(idx2,5) = round(mean(x(idx,1)),3);
            nzIdx(idx2,6) = round(mean(y(idx,1)),3);
            nzIdx(idx2,7) = round(mean(z(idx,1)),3);
        end

        % Remove points without cluster ID
        idx = nzIdx(:,2) == 0;
        nzIdx(idx,:) = [];
        
        numClust = length(sizeClust);
        combScore = zeros(numClust,1);
        allScores = zeros(1,3);
        for ii = 1:numClust
            allScores(1,1) = sizeClust(ii,1);
            allScores(1,2) = nzIdx(find(nzIdx(:,2) == ii,1),3);
            allScores(1,3) = 100-nzIdx(find(nzIdx(:,2) == ii,1),4);
            combScore(ii,1) = sum(allScores);
        end
        
        colorScore = zeros(size(nzIdx,1),1);
        for ii = 1:numClust
            colorScore(nzIdx(:,2) == ii,1) = combScore(ii,1);
        end
        
        [~,idx] = max(combScore);
        idx = find(nzIdx(:,2) == idx,1);
        optParam = nzIdx(idx,5:end);
        numObs = size(gridData.obs,1);

        subplot(3,n,i)
        scatter3(gridData.search(numObs:end,1),gridData.search(numObs:end,2),gridData.search(numObs:end,3),[],'red','filled')
        hold on
        idx = find(c >= plotThresh);
        scatter3(x(idx,1),y(idx,1),z(idx,1),[],c(idx,1),'filled')
        xlim([0,200])
        ylim([0,200])
        zlim([0,300])
        cb = colorbar;
        cb.Label.String = 'Score';
        clim(cLim);
        title(['Box ',num2str(i)])
        xlabel('Amp (uA)')
        ylabel('Feq (Hz)')
        zlabel('Dur (S)')

        subplot(3,n,i+n)
        plot(c)
        axis tight
        title('Score over time')

        subplot(3,n,i+n+2)
        scatter3(x(nzIdx(:,1),1),y(nzIdx(:,1),1),z(nzIdx(:,1),1),[],colorScore,'filled')
        xlim([0,200])
        ylim([0,200])
        zlim([0,300])
        cb = colorbar;
        cb.Label.String = 'Score';
        if isempty(clust)
            clim([0,-min(gridData.obs(:,end))]);
        else
            m = min(combScore);
            M = max(combScore);
            if m == M
                m = 0;
            end
            clim([m,M]);
        end
        title(num2str(optParam))
        xlabel('Amp (uA)')
        ylabel('Feq (Hz)')
        zlabel('Dur (S)')
    end
end