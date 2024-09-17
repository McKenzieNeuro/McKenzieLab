topDir = 'R:\ASommer\FP experiments DA-NE\Desipramine injections';
fils = getAllExtFiles(topDir,'mat',1);
kp = cellfun(@any,regexp(fils,'inject')) ;
fils = fils(kp);
[dirs,bi] = fileparts(fils);
% sl  =regexp(fils,'\');
 da  =regexp(fils,'-');
%subjs = cellfun(@(a,b) a(b(3)+1:b(4)-1),fils,sl,'uni',0);
dates = cellfun(@(a,b) a(b(1)+1:b(2)-1),fils,da,'uni',0);

%%

%define saline/desepramine sessions
k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);
DSI = [...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h4-221007-100818'} ; ...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h6-221010-120548'} ; ...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h7-221007-112247'} ; ...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h8-221014-120051'} ; ...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h9-221014-131432'} ; ...
];

SAL = [...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h4-221014-093204'};...
{'R:\ASommer\FP experiments DA-NE\Methylphenidate injections\LinearTrack_NE-220913-125322\NE2h6-221003-112857'};...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h7-221014-104711'};...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h8-221007-123704'};...
{'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h9-221007-135333'};...
];




kp_sal = ismember(dirs,SAL);
kp_DSI = ismember(dirs,DSI);
%%
    
% get all PETHs    

clear context  sessionDate subj PETH
for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        subj{i} = sessiondata.subject;
        sessionDate{i} = dates{i};
       
        ts = cell2mat(sessiondata.inject(contains(sessiondata.inject(:,1),'inject'),2));
        
        if isempty(ts)
            
            error('here')
        end
        fs =sessiondata.neural.fs_neural;
        %get PETH for each conext
        
        
        ok = nanconvn(sessiondata.neural.signal_DFoF,k');
        [ix,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts,60,3600,fs,length(sessiondata.neural.signal_DFoF));
        
         kp = ~early & ~late;
         
         ix = ix(kp,:);
         ok = ok(ix);
         ok = ok(:,1:10:end);
         PETH{i} = ok;
         ts_PETH = ts_PETH(1:10:end);
         
    end

% i
end
%subj = cellfun(@(a) a(1:5),subj,'uni',0)';
%%
% plot
PETH1  =cellfun(@(a) a(1:20:end),PETH,'uni',0);
ts_PETH1 = ts_PETH(1:20:end);

figure
plotMeanSEM(ts_PETH1,cell2mat(PETH1(kp_sal)'),'k')
hold on
plotMeanSEM(ts_PETH1,cell2mat(PETH1(kp_DSI)'),'r')
%legend('Saline','','Desipramine','')
set(gca,'fontsize',16)
%xlim([-300 3800])
ylim([-.5 7])