
fils = getAllExtFiles('R:\TransPlasticity','mat',1);
sum(contains(fils,'spikes.cellinfo.mat'))
fils = fils(contains(fils,'spikes.cellinfo.mat'));
dirs = fileparts(fils);
kp =  ~(contains(dirs,'ANI') | contains(dirs,'donot') | contains(dirs,'Thapsigargin') |  contains(dirs,'nCP1'));
dirs = dirs(kp);
%%

for i = 1:length(dirs)
    cd(dirs{i})
    dirs{i}
    fils = getAllExtFiles( dirs{i},'mat',0);
    if ~any(contains(fils,'.cell_metrics.cellinfo.mat'))
        
        session = sessionTemplate(pwd, 'showGUI', false);   % or load your session struct
        cell_metrics = ProcessCellMetrics('session', session,'showGUI', false,'manualAdjustMonoSyn',false,'showFigures',false,...
            'excludeMetrics',{'monoSynaptic_connections'});
    end
close all    
    
end