
topDir = 'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field';
fils = getAllExtFiles(topDir,'mat',1);
kp = cellfun(@any,regexp(fils,'inject')) ;
fils = fils(kp);
[dirs,bi] = fileparts(fils);
% sl  =regexp(fils,'\');
 da  =regexp(fils,'-');
% subjs = cellfun(@(a,b) a(b(3)+1:b(4)-1),fils,sl,'uni',0);
dates = cellfun(@(a,b) a(b(1)+1:b(2)-1),fils,da,'uni',0);
%%


% save sessiondata struct
for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
   % if ~exist('sessiondata.mat')
  
        
        saveAllNeuralPosition(pwd)
        
  %  end
    
      v=load('inject.mat');
    load('sessiondata.mat')
    sessiondata.inject =v.data;
   
  
    save('sessiondata','sessiondata','-v7.3')
    i
end



%%

%define saline/yohimbine sessions
k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);
saline = [...
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221004-105054'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h7-221004-114833'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h6-221011-125645'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h9-221011-144957'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221011-135518'}
    ];

YOH = [...
%    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221005-102443'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h7-221005-111958'}
    %{'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221011-135518'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h6-221012-091931'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221012-101436'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h9-221012-122213'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221017-150806'}
    ];


kp_sal = ismember(dirs,saline);
kp_YOH = ismember(dirs,YOH);
%%
    
% get all PETHs    
kp_sal = ismember(dirs,saline);

clear context  sessionDate subj PETH
for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        subj{i} = sessiondata.subject;
        sessionDate{i} = dates{i};
       
        ts = cell2mat(sessiondata.inject(contains(sessiondata.inject(:,1),'novel'),2));
        
        if isempty(ts)
            
            error('here')
        end
        fs =sessiondata.neural.fs_neural;
        %get PETH for each conext
        
        
        ok = nanconvn(sessiondata.neural.signal_DFoF,k');
        ts_neural = (1:length(ok))/fs;
       % [~,ixx2] = histc(sessiondata.inject{1,2}-100,ts_neural);
       % [~,ixx1] = histc(sessiondata.inject{1,2}-600,ts_neural);
       % baseline = nanmean(ok(ixx1:ixx2));
      %  ok = ok-baseline;
        [ix,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts,60,60,fs,length(sessiondata.neural.signal_DFoF));
        
         kp = ~early & ~late;
         
         ix = ix(kp,:);
         
         PETH{i} = ok(ix);
         
    end

i
end
subj = cellfun(@(a) a(1:5),subj,'uni',0)';
%%
PETH1  =cellfun(@(a) a(1:10:end),PETH,'uni',0);
ts_PETH1 = ts_PETH(1:10:end);
% plot
close all

figure
plotMeanSEM(ts_PETH1,cell2mat(PETH1(kp_sal)'),'r')
ylim([-1 2])
f = gcf ;
f.Renderer = 'painters'
%saveas(f,'novel_sal','epsc');
%figure
plotMeanSEM(ts_PETH1,cell2mat(PETH1(kp_YOH)')-.5,'k')
ylim([-1 2])
f = gcf ;
f.Renderer = 'painters'
%saveas(f,'novel_yoh','epsc');