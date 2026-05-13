%Projective Transformation

movingPoints = [edg(:,1), edg(:,2)];
if contains(sessiondata.contextEntry{j,1}, 'home' )
    fixedPoints = [0 0; 24.88 0; 24.88 14.85; 0 14.85];
elseif contains(sessiondata.contextEntry{j,1}, 'ontext1')
    fixedPoints = [0 0; 24.88 0; 24.88 14.85; 0 14.85];
elseif contains(sessiondata.contextEntry{j,1}, 'ontext2')
    fixedPoints = [0 0; 39.43 0; 39.43 18.1; 0 18.1];
elseif contains(sessiondata.contextEntry{j,1}, 'BB')
    fixedPoints = [0 0; 16.47 0; 16.47 16.47; 0 16.47];
elseif contains(sessiondata.contextEntry{j,1}, 'CB')
    fixedPoints = [0 0; 9.97 0; 9.97 9.97; 0 9.97];
elseif contains(sessiondata.contextEntry{j,1}, 'RB')
    fixedPoints = [0 0; 14.67 0; 14.67 9.90; 0 9.90];
elseif contains(sessiondata.contextEntry{j,1}, 'CC')
    fixedPoints = [0 0; 13.77 0; 13.77 13.77; 0 13.77];
elseif contains(sessiondata.contextEntry{j,1}, 'OTT')
    fixedPoints = [0 0; 25.60 0; 25.60 15.60; 0 15.60];
elseif contains(sessiondata.contextEntry{j,1}, 'WB')
    fixedPoints = [0 0; 10.00 0; 10.00 10.00; 0 10.00];
end
tform = fitgeotrans(movingPoints , fixedPoints , "projective");

[edg(:,1),edg(:,2)] = transformPointsForward(tform, edg(:,1)',edg(:,2)');
[LE(:,1),LE(:,2)] = transformPointsForward(tform, LE(:,1)', LE(:,2)');
[RE(:,1),RE(:,2)] = transformPointsForward(tform, RE(:,1)', RE(:,2)');
[TB(:,1),TB(:,2)] = transformPointsForward(tform, TB(:,1)', TB(:,2)');
