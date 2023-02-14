
dirs{1} = 'R:\McKenzieLab\DANEHippocampalResponse\DA2h2\Novel Environment\DA2h2-210817-141645';
dirs{2} = 'R:\McKenzieLab\DANEHippocampalResponse\DA2h2\Novel Environment\DA2h2-210818-151936';

neural =[];
vel =[];
acc = [];
kp_Novel =[];
sesID = [];
for i = 1:2

    cd(dirs{i})


    %load struct with neural/behavioral data
    load('sessiondata.mat')

    %get context transitions
    data = sessiondata.contextEntry;

    %get times for each context entry
    context_entry = cell2mat(data(:,2));

    %sort and get which index goes in which order (b)
    [context_entry,b] = sort(context_entry);

    %sort the edge/transition by time
    data = data(b,:);

    %define entry times
    epochs_on  = context_entry;

    %define exit times
    epochs_off = [context_entry(2:end); sessiondata.behavior.ts_video(end)];

    %build matrix of onsets and offsets
    epochs = [epochs_on epochs_off];

    %find epochs in home cage
    home = cellfun(@any,regexp(data(:,1),'home'));
    inNovel = ~home;


    %get times in home cage and in novel context
    ts_video = sessiondata.behavior.ts_video;

    kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
    kp_homet = cell2mat(kp_homet');
    kp_homet = any(kp_homet,2);

    kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
    kp_Novelt = cell2mat(kp_Novelt');
    kp_Novelt = any(kp_Novelt,2);
    
    %get velocity
    velt = sessiondata.behavior.vel;
    acct = sessiondata.behavior.acc;

    %downsample neural data to video clock
    neuralt = sessiondata.neural.signal_DFoF;
    fs = sessiondata.neural.fs_neural;

    ts_neural = (1:length(neuralt))/fs;

    neuralt = double(interp1(ts_neural,neuralt,ts_video));

    nSamples = length(neuralt);
    % concatenate each session
    neural = [neural;neuralt];
    acc = [acc;acct];
    vel = [vel;velt];
    kp_Novel = [kp_Novel;kp_Novelt];
    sesID = [sesID; i *ones(nSamples,1)];

i
end
%%

dat = table(neural,vel,acc,double(kp_Novel),sesID,'VariableNames',{'neural','velocity','acceleration','kp_novel','sesID'});


glme = fitglme(dat,'neural ~ 1 + kp_novel + velocity + acceleration + (1|sesID)');
