fils = getAllExtFiles(pwd,'mat',1);

%find all linear track files
kp = cellfun(@any,regexpi(fils,'trackON'));
%nokp = ~cellfun(@any,regexpi(fils,'SOR'));
fils = fils(kp);

kp = false(length(fils),1);
for i = 1:length(fils)
    
    [a,b] = fileparts(fils{i});
    cd(a)
    
    data = TDTbin2mat(a);
    
    if isfield(data.streams,'x500D')
        
        streams{1} = 'x500D';
        streams{2} = 'x450D';
    elseif  isfield(data.streams,'x465A')
        streams{1} = 'x465A';
        streams{2} = 'x405A';
    end
    fs = data.streams.(streams{1}).fs;
    
    
    ts = (1:length(data.streams.(streams{1}).data))/fs;
    plot(ts,data.streams.(streams{1}).data,'k')
    hold on
    plot(ts,data.streams.(streams{2}).data,'r')
    
    
    kp(i) = input('1=good; 0 = bad ');
    close all
    
end

%%

subjID = cellfun(@(a,b) a(b(2)+1:b(3)-1), fils,regexp(fils,'\'),'uni',0);
sesDate = cellfun(@(a,b) a(b(1)+1:b(2)-1), fils,regexp(fils,'-'),'uni',0);
dirN = fileparts(fils);


%%
%video code good files



[a,b] = fileparts(fils{i}); cd(a)
sm_labelTDTvideo
%%





[data,ix,signal_smoothed] = sm_QC_DANE(dirName,sensor)
