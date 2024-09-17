

%define saline/yohimbine sessions
k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);
dirs = [...
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221004-105054'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h7-221004-114833'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h6-221011-125645'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h9-221011-144957'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221011-135518'}
    
    %    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221005-102443'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h7-221005-111958'}
    %{'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221011-135518'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h6-221012-091931'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h8-221012-101436'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h9-221012-122213'}
    {'R:\ASommer\FP experiments DA-NE\Yohimbine - Novel Field\NE2h4-221017-150806'}
    
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h8-221007-123704'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h9-221007-135333'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h4-221014-093204'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h7-221014-104711'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h4-221007-100818'} ;...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h6-221010-120548'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h7-221007-112247'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h8-221014-120051'};...
    {'R:\ASommer\FP experiments DA-NE\Desipramine injections\NE-Desipramine_injections-221007-100452\NE2h9-221014-131432'};...
    ];

k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);


%%



for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        subj{i} = sessiondata.subject;

       
        ts = cell2mat(sessiondata.inject(contains(sessiondata.inject(:,1),'inject'),2));
        
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
        [ix,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts,300,300,fs,length(sessiondata.neural.signal_DFoF));
        
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