masterDir = 'R:\McKenzieLab\DANEHippocampalResponse';
fils = getAllExtFiles(masterDir,'mat',1);

%kp_edges = cellfun(@any,regexp(fils,'arena_edges'));
%kp_trans = cellfun(@any,regexp(fils,'transition'));
kp_Novel = cellfun(@any,regexp(fils,'Novel Env'));

fils_N = fils(kp_Novel);
[dirs] = fileparts(fils_N);

dirs = unique(dirs);


count = 0;

%%
for i = 1:length(dirs)
    i
    % change directory into each (ith) directory in list
    cd(dirs{i});
    
    %check if we have tracking/neural data
    if exist('sessiondata.mat')
        
        %load tracking data
        load('sessiondata.mat');
        
        %load context edges
        load('contextTransition1.mat');
        
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
        
        
        %plot
        figure
        
        %define number of plots needed for even numbers of rows/cols
        N = ceil(sqrt(size(epochs,1)));
        ax  = tight_subplot(N,N);
        
        %loop over contexts
        for j = 1:size(epochs,1)
            
            %define which axis to plot
            axes(ax(j))
            
            %select the subset of times (in the video) between the onset
            %and offset of the jth epoch
            kp = sessiondata.behavior.ts_video > epochs(j,1) & sessiondata.behavior.ts_video<epochs(j,2);
            
            %plot tracking data for those subsets of time (kp)
            plot(sessiondata.behavior.position.left_ear(kp,1),sessiondata.behavior.position.left_ear(kp,2),'.');
            
            
            %don't overwrite your plot!
            hold on
         
            %plot the context edges (flipped around the vertical dimension)
            plot(data{j,3}(:,1),480-data{j,3}(:,2), 'x')
            
        end
        
        %wait for the user to press a button
        waitforbuttonpress
        
        %clear the plot
        close all
    else
        
        %signal (for Infania) that we need to get the tracking data for these sessions
        disp([dirs{i} ' not loaded into sessiondata'])
    end
end

