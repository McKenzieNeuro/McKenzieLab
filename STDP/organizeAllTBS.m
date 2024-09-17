% define directories to analyze
masterDir = 'E:\Dropbox\UNM\data\';

dirs = [ ...
    
{[ masterDir 'STDP\STDP4\STDP4_210625']} ; ...%loaded

{[ masterDir 'STDP\STDP4\STDP4_210627']} ; ...%loaded
{[ masterDir 'STDP\STDP4\STDP4_210628']} ; ...%loaded
{[ masterDir 'STDP\STDP4\STDP4_210701']} ; ...%loaded
{[ masterDir 'STDP\STDP4\STDP4_210702']} ; ...%loaded
{[ masterDir 'STDP\STDP4\STDP4_210703']} ; ...%loaded
{[ masterDir 'STDP\STDP4\STDP4_210706']} ; ...
{[ masterDir 'STDP\STDP4\STDP4_210707']} ; ...
{[ masterDir 'STDP\STDP4\STDP4_210708']} ; ...

];


%%
% initialize all variables
DR_post_all =[];
DR_pre_all =[];
isNonStim_all =[];
isPaired_all =[];
isunPaired_all =[];
shankdiff_all =[];


basename = 'amplifier_analogin_auxiliary_int16';


%loop through directories
for i = 1:length(dirs)
    
    
    %go to directory and load all files
    cd(dirs{i})
    %basename = bz_BasenameFromBasepath(dirs{i});
    
    load([basename '.spikes.cellinfo.mat'])
    load('pulse_info.mat')
    
    
    %load([basename '.CellParams.mat'])
    %load([basename '-states.mat'])
    %load('mono_res.mat')
    
    % define beginning and end of each DR
    pre1_start = pulse_info(1).doseResponse(1,1);
    pre1_end = pulse_info(1).doseResponse(1,2);
    post1_start = pulse_info(1).doseResponse(2,1);
    post1_end = pulse_info(1).doseResponse(2,2);
    
    
    %loop through shanks
    for j =1:4
        % get pre1 pulses
        pre1 = InIntervals(pulse_info(j).time(:,1),[pre1_start pre1_end]);
        post1 = InIntervals(pulse_info(j).time(:,1),[post1_start post1_end]);
        
        
        
        pulse{j} = [pulse_info(j).time(:,1) pulse_info(j).time(:,2) pre1 post1 pulse_info(j).id];
        
    end
    
    %build dose response for each blue light
    
    % find which is blue
    
    kp = find(contains({pulse_info.color},'B'));
    
    %loop over blue shanks
    clear DR_pre DR_post
    for j = kp
        pre_ts =  pulse{j}( pulse{j}(:,3)==1,1:2);
        pre_id = pulse{j}( pulse{j}(:,3)==1,5);
        post_ts =  pulse{j}( pulse{j}(:,4)==1,1:2);
        post_id =  pulse{j}( pulse{j}(:,4)==1,5);
        
        %assume DR are 100ms (for now)
        [binnedPop] = populationMatrix(spikes,0.002,.095,1,pre_ts(:,1));
        
        
        for k = 1:size(binnedPop,1)
            DR_pre{j}(k,:) = accumarray(pre_id,squeeze(binnedPop(k,:,:)),[ 5 1],@nanmean,nan);
        end
        [binnedPop] = populationMatrix(spikes,0.002,.095,1,post_ts(:,1));
        
        
        for k = 1:size(binnedPop,1)
            DR_post{j}(k,:) = accumarray(post_id,squeeze(binnedPop(k,:,:)),[ 5 1],@nanmean,nan);
        end
        
        
        
        
    end
    
    %find red shanks
    kp = find(contains({pulse_info.color},'R'));
    paired_dac = find(pulse_info(1).pairing(kp,:)==1);
    paired_dac = paired_dac(:)'; % force to be a row
    paired_shank =[];
    for d = paired_dac
        
        paired_shank =  [paired_shank;pulse_info(d).shank];
        
    end
    
    isPaired = ismember(spikes.shankID,paired_shank);
    
    unpaired_dac = find(pulse_info(1).pairing(kp,:)==-1);
    unpaired_dac = unpaired_dac(:)'; % force to be a row
    unpaired_shank =[];
    for d = unpaired_dac
        
        unpaired_shank =  [unpaired_shank;pulse_info(d).shank];
        
    end
    
    isunPaired = ismember(spikes.shankID,unpaired_shank);
    
    if any(unpaired_shank) | any(paired_shank)
        nonstim  = setdiff(1:4,[unpaired_shank;paired_shank]);
        
    else
        %control day without any red light
    end
    isNonStim = ismember(spikes.shankID,nonstim)';
    
    % loop over DR shanks and concat into matrix with distance from
    % recording and stim shank
    
    for j = 1:length(DR_pre)
        
        if ~isempty(DR_pre{j})
            
            %calculate the distance between jth stim shank and neuron
            shankdiff = abs(pulse_info(j).shank-spikes.shankID)';
            DR_pretmp = DR_pre{j};
            DR_posttmp = DR_post{j};
            
            DR_pre_all = [DR_pre_all;DR_pretmp];
            DR_post_all = [DR_post_all;DR_posttmp];
            shankdiff_all = [shankdiff_all;shankdiff(:)];
            isunPaired_all = [isunPaired_all;isunPaired(:)];
            isPaired_all = [isPaired_all;isPaired(:)];
            isNonStim_all = [isNonStim_all;isNonStim(:)];
            
        end
    end
end


%%
figure
% plot

% change in DR on same shank for paired
plot(nanmean(DR_post_all(shankdiff_all==0 & isPaired_all==1,:) - DR_pre_all(shankdiff_all==0 & isPaired_all==1,:)))
hold on

% change in DR on same shank for unpaired
plot(nanmean(DR_post_all(shankdiff_all==0 & isunPaired_all==1,:) - DR_pre_all(shankdiff_all==0 & isunPaired_all==1,:)))

