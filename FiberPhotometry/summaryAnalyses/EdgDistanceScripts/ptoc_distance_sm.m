
function distanceCircle = ptoc_distance_sm(edges, points)

npoints = size(points, 2);
distanceCircle = nan(npoints, 1);

centX = range(edges(:,1))/2 + min(edges(:,1));
centY = range(edges(:,2))/2 + min(edges(:,2));
dCenEdg = mean(range(edges)/2 );

X = points(:,1);
Y = points(:,2);

X = X  - centX;
Y = Y -  centY;
dCenLear = sqrt(X.^2 + Y.^2);
for j = 1:length(X)
    
    %%finds distance from center pts to body coord
    % cent2lear = [X,centX ; Y,centY];
    
    %%finds distance from body coord. to edge
    
    distanceCircle(j)  = abs(dCenLear(j) - dCenEdg);
    %        distanceCircle(j,1) = dLearEdg(j,1);
    
end

end