function epochs = excludeEpochs(goodEpoch, badEpoch)

% epochs = Utils.ExcludeEpochs(goodEpochs, badEpochs)

%

% Removes periods of time from goodEpochs according to badEpochs

%

% Contributed by Sam Mckenzie on 12/13/2011

% Comments and addition to CMBH Toolbox: 12/16/2011 arb

 

[~, order] = sort(goodEpoch(:,1), 'ascend'); % sort good and bad epochs by start time

 

goodEpoch = goodEpoch(order, :);

 

[~, order] = sort(badEpoch(:,1), 'ascend');

 

badEpoch = badEpoch(order, :);

 

newEpoch=[];

for i=1:size(badEpoch,1)

    change1= badEpoch(i,1) > goodEpoch(:,1) & badEpoch(i,1)<goodEpoch(:,2);

    newEpoch=[newEpoch; goodEpoch(change1,1) badEpoch(i,1).*ones(sum(change1),1)];

    goodEpoch(change1,:)=[badEpoch(i,2).*ones(sum(change1),1)  goodEpoch(change1,2)];

    goodEpoch=goodEpoch(goodEpoch(:,1)<goodEpoch(:,2),:);

    change2= badEpoch(i,2) > goodEpoch(:,1) & badEpoch(i,2)<goodEpoch(:,2);

    goodEpoch(change2,:)=[badEpoch(i,2).*ones(sum(change2),1)  goodEpoch(change2,2)];

    change3=badEpoch(i,1)<goodEpoch(:,1) & badEpoch (i,2) > goodEpoch(:,2);

    goodEpoch=goodEpoch(~change3,:);

end

epochs=[goodEpoch;newEpoch];

epochs=sortrows(epochs);

 

end