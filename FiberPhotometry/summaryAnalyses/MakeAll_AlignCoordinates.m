fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

[dirs] = fileparts(fils);

dirs =  unique(dirs);
%%
for i = 1:length(dirs)
    cd(dirs{i})
    i
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


        nSamples = length(sessiondata.behavior.position.left_ear);
        sessiondata.behavior.position.left_ear_cor = nan(nSamples,2);
        sessiondata.behavior.position.right_ear_cor = nan(nSamples,2);
        sessiondata.behavior.position.tail_base_cor = nan(nSamples,2);
        sessiondata.behavior.position.context = cell(nSamples,1);
        %loop over contexts
        for j = 1:size(epochs,1)



            %select the subset of times (in the video) between the onset
            %and offset of the jth epoch
            kp = sessiondata.behavior.ts_video >= epochs(j,1) & sessiondata.behavior.ts_video<epochs(j,2);

            %keep data for those subsets of time (kp)
            LE =  sessiondata.behavior.position.left_ear(kp,:);
            RE =  sessiondata.behavior.position.right_ear(kp,:);
            TB = sessiondata.behavior.position.tail_base(kp,:);

            %define the edges of the context for that subst of time (note
            %we need to flip by 480 since the edges are defined in flipped
            %coordinates  (for some reason) and the height of the video is 480px
            
            edg = [data{j,3}(:,1) 480-data{j,3}(:,2)];
            
            %find the left edge
            [~,b] = sort(edg(:,1));
            left = edg(b(1:2),:);
            
            %find the bottom left point
            [~,b] = min(left(:,2));
            BL = left(b,:);

            % translate the position and edges so the bottom left is at the
            % origin


            LE(:,1) = LE(:,1) - BL(1);
            LE(:,2) = LE(:,2) - BL(2);

            RE(:,1) = RE(:,1) - BL(1);
            RE(:,2) = RE(:,2) - BL(2);

            TB(:,1) = TB(:,1) - BL(1);
            TB(:,2) = TB(:,2) - BL(2);

            edg(:,1) = edg(:,1) - BL(1);
            edg(:,2) = edg(:,2) - BL(2);

            
            %now find top left
            [~,b] = sort(edg(:,1));
            left = edg(b(1:2),:);
            [~,b] = max(left(:,2));
            TL = left(b,:);
            
            % get angle of offset
            deg_offset = (atan(TL(1)/TL(2)));
            
            %find the rotation matrix to offset that angle
            R = rot2d( -deg_offset );
         
            %rotate the position data
            LE = LE*R;
            RE = RE*R;
            TB = TB*R;


            %rotate the edges
            ok =(edg*R);

            
            
            %keep long edge up/down
            xyrange = range(ok);

            if xyrange(1)>xyrange(2)
                R = rot2d( 90 );
               
                LE = LE*R;
                RE = RE*R;
                TB = TB*R;
                ok =(ok*R);
                
                %shift over by new width
                w = range(ok(:,1));
                
                LE(:,1) = LE(:,1) + w;
                RE(:,1) = RE(:,1) + w;
                TB(:,1) = TB(:,1) + w;
                ok(:,1) = ok(:,1) +w;
                
            end
            
            % rescale by the long edge (y-axis)
            LE = LE/range(ok(:,2));
            RE = RE/range(ok(:,2));
            TB = TB/range(ok(:,2));

            sessiondata.behavior.position.left_ear_cor(kp,:) = LE;
            sessiondata.behavior.position.right_ear_cor(kp,:) = RE;
            sessiondata.behavior.position.tail_base_cor(kp,:) = TB;
            sessiondata.behavior.position.context(kp) = data(j,1);




           
        end
        sessiondata.contextEntry = data;
        
        save('sessiondata.mat','sessiondata','-v7.3')

    end
end