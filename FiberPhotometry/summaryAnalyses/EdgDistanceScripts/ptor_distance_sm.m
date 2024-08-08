%function distanceRectangle = ptor_distance_sm(edges, points)
%%%FOR HOME CONTEXT OR RECTANGULAR %%
%%%plot the perimeter of the context
%npoints = size(points, 1);
%distEdgAll = zeros(npoints,4);
%smalDist = zeros(npoints,1);

%x1 = [edges(:,1)];
%x1(end+1) = edges(1:1);
%y1 = [edges(:,2)];
%y1(end+1) = [edges(1,2)];
%distanceRectangle = nan(size(points,1),1);
%%find distance from the context (assuming edges are aligned
%% each edge shares either a common x-axis or y-axis)
%for m = 1:length(points)
%    dist2side = min(edges(:,1) -  points(1,1));
%    dist2topBottom = min(edges(:,2) -  points(1,2));
%    
%    distanceRectangle(m,1) = min(dist2side,dist2topBottom);
%end
%
%end


function distanceRectangle = ptor_distance_sm(edges, points)
%%FOR HOME CONTEXT OR RECTANGULAR %%
%%plot the perimeter of the context
npoints = size(points, 1);


x1 = [edges(:,1)];
x1(end+1) = edges(1:1);
y1 = [edges(:,2)];
y1(end+1) = [edges(1,2)];
distanceRectangle = nan(size(points,1),1);
%find distance from the context (assuming edges are aligned
% each edge shares either a common x-axis or y-axis)
    for m = 1:npoints
       
       dist2side = min(abs(edges(:,1) -  points(m,1)));
       dist2topBottom = min(abs(edges(:,2) -  points(m,2)));
   
       distanceRectangle(m,1) = min(dist2side,dist2topBottom);
    end

end