% This function calculates the distance between a set of points 
% and an ellipse defined by a set of edges. The ellipse is defined 
% by its center point, half-width and half-height. 
% The function first calculates the center and dimensions of the ellipse 
% using the edges input, and then plots the ellipse using the x and y 
% coordinates of the points on the ellipse's perimeter. It then loops 
% through each point in the points input and calculates the Euclidean 
% distance between that point and every point on the ellipse's perimeter. 
% The minimum distance is then stored in the distanceEllipse output 
% for that point.

function distanceEllipse = ptoe_distance_sm_progress(edges, points)
    npoints = size(points, 1);
    distanceEllipse = zeros(npoints, 1);
    eclipse_center_x = 0.5 * (max(edges(:, 1)) + min(edges(:, 1)));  % Ellipse center in x
    eclipse_center_y = 0.5 * (max(edges(:, 2)) + min(edges(:, 2)));  % Ellipse center in y
    half_width = 0.5 * (max(edges(:, 1)) - min(edges(:, 1)));  % Ellipse half-width
    half_height = 0.5 * (max(edges(:, 2)) - min(edges(:, 2)));  % Ellipse half-height
    
    a = half_width; % horizontal radius
    b = half_height; % vertical radius
    x0 = eclipse_center_x; % x0,y0 ellipse centre coordinates
    y0 = eclipse_center_y;
    t = -pi:0.01:pi;
    x = x0+a*cos(t);
    y = y0+b*sin(t);
    plot(x,y)
    
    numEllipseP = length(x);
    distPoint = zeros(npoints,1);
    for i = 1:npoints
        
        distEllipse = (zeros(numEllipseP,1));
        for ii = 1:numEllipseP
            
            distEllipse(ii,1) = pdist([points(i,1),points(i,2);x(1,ii),y(1,ii)],'euclidean');
        end
        distanceEllipse(i, 1) =  min(distEllipse);
    end
%     
%     numEllipseP = length(x);
% %     npoints = size(points,1);
% %     distanceEllipse = zeros(npoints,1);
%     numQuad = 3;
%     numIter = 5;
% 
%     for i = 1:npoints
%         
%         start = 1;
%         numP = numEllipseP;
%         for ii = 1:numIter
%             
%             step = round(numP/numQuad);
%             distQuad = zeros(numQuad+1,2);
%             distQuad(1,1) = start; 
%             for iii = 2:numQuad
%                 distQuad(iii,1) = start+step*(iii-1);        
%             end
%             distQuad(end,1) = numEllipseP; 
% 
%             for iii = 1:numQuad+1    
%                 distQuad(iii,2) = pdist([points(i,1),points(i,2);x(1,distQuad(iii,1)),y(1,distQuad(iii,1))],'euclidean');
%             end
% 
%             [~,iDist] = sort(distQuad(:,2));
% 
%             Q(1,1) = distQuad(iDist(1,1),1);
%             Q(2,1) = distQuad(iDist(2,1),1);
%             Q = sort(Q);
% 
%             if Q(1,1) == 1 && Q(2,1) == numEllipseP
%                 break
%             end
% 
%             start = Q(1,1);
%             numP = length(Q(1,1):Q(2,1));
%         end
% 
%         if Q(1,1) == 1 && Q(2,1) == numEllipseP
%             distanceEllipse(i,1) = min(distQuad(:,2));
%         else
%             distEllipse = zeros(numP,1);
%             count = 1;
%             for ii = Q(1,1):Q(2,1)
%                 distEllipse(count,1) = pdist([points(i,1),points(i,2);x(1,ii),y(1,ii)],'euclidean');
%                 count = count+1;
%             end    
%             distanceEllipse(i,1) =  min(distEllipse);
%         end
%     end    
%      
% %      binX = 0:.05:2; % this is the resolution to bin the rows
% %      binY = -1.5:0.05:0.5; % this is the resolution to bin the columns
% %      distan = totalEdgDist(kp);
% %      [n,~,~,ix] = histcn(points,binX,binY); % bin your 2D data
% %      kp  = all(ix>0,2); % exclude all points that are outside of your bins
% %      mean_dist = accumarray(ix(kp,:),distPoint(kp,:),[],@nanmean,nan); % calculate the mean distance from the edge at every binned position
% %      figure
% %      imagesc(mean_dist)
% 
end


