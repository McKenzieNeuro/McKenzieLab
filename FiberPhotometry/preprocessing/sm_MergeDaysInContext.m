function sm_MergeDaysInContext(topDir)
%
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'sessiondata.mat');
fils = fils(kp);
%%
%get data and sort
clear date
for i = 1:length(fils)
    dirN = fileparts(fils{i});
    [~,subj] = fileparts(dirN);
    subj_split = split(subj, "-");
    
    datet = subj_split{2};
    %datet = [datet(3:4) '/' datet(1:2) '/'  datet(5:6)];
    date(i) = str2num(datet);
    
    
end

[~,ord] = sort(date);

fils = fils(ord);
%%
%get all contexts
for i =1 :length(fils)
    load(fils{i})
    if i ==1
        
        hc = contains(sessiondata.contextEntry(:,1),'home');
        sessiondata.contextEntry(~hc,5) = {0};
        allcon = unique(sessiondata.contextEntry(~hc,1));
        numvis = zeros(length(allcon),1);
    else
        hc = contains(sessiondata.contextEntry(:,1),'home');
        uCon = unique(sessiondata.contextEntry(~hc,1));
        
        for j = 1:length(uCon)
            [ix,b] = ismember(uCon(j),allcon);
            idx = contains(sessiondata.contextEntry(:,1),uCon);
            if ix==1
                numvis(b) = numvis(b)+1;
                sessiondata.contextEntry(idx,5) = {numvis(b)};
            else
                allcon = [allcon;uCon(j)];
                numvis = [numvis;0];
                 sessiondata.contextEntry(idx,5) = {0};
            end
        end
    end
    save(fils{i},'sessiondata')
end
%%
end