ix = 1;
for i = 1:length(dirs)
    cd(dirs{i})
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        N(ix) =   min(length(sessiondata.behavior.rewardTimeLeft),length(sessiondata.behavior.rewardTimeRight));
        dur(ix) =   range([sessiondata.behavior.rewardTimeLeft;sessiondata.behavior.rewardTimeRight]);
        ix = ix+1
    end
end