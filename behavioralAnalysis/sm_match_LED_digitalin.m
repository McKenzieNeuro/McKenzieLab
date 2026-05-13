function sm_match_LED_digitalin(outdir,up,dwn,whls,dur,fs,fils,varargin)

p = inputParser;
addParameter(p,'merge_cond',[],@isstr);
addParameter(p,'threshF',[],@isnumeric);

parse(p,varargin{:});

merge_cond = p.Results.merge_cond;
threshF =  p.Results.threshF;

%max time between blinks (s)
LED_ISI_max = 10;


% clean up digitalin

kp = dwn-up>.15; % pulses must be > 150ms
up = up(kp);
dwn = dwn(kp);


if ~iscell(whls)
    whls{1} = whls;
end


if ~isempty(merge_cond) && ~strcmp(merge_cond,'1Intan1Vid')
    
    switch merge_cond
        case '2Intan2Vid'
            
            if length(whls)>1
                
                %find subset of ups for first video
                edg = cumsum([0 dur]);
            end
            
        case '1Intan2Vid'
            
            %here we will need to infer the start of the new video by a
            %large gap
            idx = find(diff(up)>LED_ISI_max);
            edg = cumsum([0 up(idx)+.1 max(up)+10]);
            
    end
    [~,vid_num] = histc(up,edg);
    
else
    vid_num = ones(length(up),1);
    
end

numVids = max(vid_num);


if  numVids ~= length(whls)
    error('cannot match pulses with videos')
end



for v = 1:numVids
    
    %get the LEDs and pulses for each video
    
    LED = whls{v};
    dwns = dwn(vid_num==v);
    ups = up(vid_num==v);
    if  isempty(threshF)
        threshF = prctile(diff(LED),98);
    end
    
    
    up_vid = find(diff([0;[0;diff(LED)>threshF]])>0);
    dwn_vid = find(diff([0;[0;diff(LED)<-threshF]])>0);
    
    dt = median(dwns-ups);
    dt_fr = median(round(dt*fs));
    
    kp = abs((up_vid+dt_fr)-bestmatch(up_vid+dt_fr,dwn_vid))<5;
    up_vid = up_vid(kp);
    
    
    
    max_dt  =max(diff(ups))+.2;
    %%
    % now find matches
    minJ = 1;ts_syn =[];
    it=1;
    idx=1;
    
    allfound = false;
    
    while ~allfound
        
        if it<length(ups)-3 && minJ+3 <= length(up_vid)
            dt3 = diff(ups(it:(it+3)));
            dt3_vid = diff(up_vid(minJ:minJ+3))/fs;
            
            if all(abs(dt3 - dt3_vid)<.15)
                
                ts_syn = [ts_syn;ups(it) up_vid(minJ) up_vid(minJ)/fs minJ];
                
                
                minJ = minJ+1;
                it=it+1;
                
            else
                
                
                %find next set of videos that match
                %loop through all videos and find next set of digitalin
                %that matched
                
                found = false;
                
                if ~isempty(ts_syn)
                    minJ = ts_syn(end,4);
                else
                    minJ=1;
                end
                for jj = minJ+1:length(up_vid)-4
                    
                    if ~found
                        dt3_vid = diff(up_vid(jj:jj+4))/fs;
                        for ii = it:length(ups)-4
                            dt3 = diff(ups(ii:(ii+4)));
                            if all(abs(dt3 - dt3_vid)<.15)
                                
                                % check if it is on regression line
                                
                                if size(ts_syn,1)>=3
                                    kp=~isnan(ts_syn(:,2));
                                    ok = fit(ts_syn(kp,1) , ts_syn(kp,2),'poly1');
                                    pt =feval(ok,ups(ii));
                                end
                                if size(ts_syn,1)<3 || abs(up_vid(jj)-pt)<500
                                    it = ii;
                                    minJ = jj;
                                    found = true;
                                    break
                                else
                                    er(idx) = abs(up_vid(jj)-pt);
                                    idx = idx+1;
                                end
                            end
                        end
                        
                        if ii==length(ups)-4 && jj ==length(up_vid)-4
                            it = ii+1;
                            break
                        end
                        
                    elseif found
                        break
                        
                        
                    end
                    
                end
            end
        else
            allfound = true;
        end
    end
    
    
    kp = abs(diff(ts_syn(:,1))-diff(ts_syn(:,3)))<1;
    ts_syn = ts_syn(kp,:);
    %%
    % try reverse
    if size(ups,1) ~= size(ts_syn,1)
        
        
        
        % now find matches
        minJ = length(up_vid)-3;
        it=length(ups)-3;
        idx=1;
        ts_syn_back =[];
        allfound = false;
        
        
        while ~allfound
            
            if it>0
                dt3 = diff(ups(it:(it+3)));
                dt3_vid = diff(up_vid(minJ:minJ+3))/fs;
                
                if all(abs(dt3 - dt3_vid)<.15)
                    
                    ts_syn_back = [ts_syn_back;ups(it) up_vid(minJ) up_vid(minJ)/fs minJ];
                    
                    
                    minJ = minJ-1;
                    it=it-1;
                    
                else
                    
                    
                    %find next set of videos that match
                    %loop through all videos and find next set of digitalin
                    %that matched
                    
                    found = false;
                    
                    if ~isempty(ts_syn_back)
                        minJ = ts_syn_back(end,4);
                    else
                        minJ=length(up_vid)-4;
                    end
                    if minJ>length(up_vid)-4
                        minJ =  length(up_vid)-4;
                    end
                    for jj = minJ+1:-1:1
                        
                        if ~found
                            dt3_vid = diff(up_vid(jj:jj+3))/fs;
                            for ii = it:-1:1
                                dt3 = diff(ups(ii:(ii+3)));
                                if all(abs(dt3 - dt3_vid)<.15)
                                    
                                    % check if it is on regression line
                                    
                                    if size(ts_syn_back,1)>=3
                                        kp=~isnan(ts_syn_back(:,2));
                                        ok = fit(ts_syn_back(kp,1) , ts_syn_back(kp,2),'poly1');
                                        pt =feval(ok,ups(ii));
                                    end
                                    if size(ts_syn_back,1)<3 || abs(up_vid(jj)-pt)<500
                                        it = ii;
                                        minJ = jj;
                                        found = true;
                                        break
                                    else
                                        er(idx) = abs(up_vid(jj)-pt);
                                        idx = idx+1;
                                    end
                                end
                            end
                            
                            if ii==1 && jj ==1
                                it = 0;
                                break
                            end
                            
                        elseif found
                            break
                            
                            
                        end
                        
                    end
                end
            else
                allfound = true;
            end
        end
        
        
        
        
    end
    kp = (~ismember(ts_syn_back(:,1),ts_syn(:,1)));
    ts_syn_back = ts_syn_back(kp,:);
    [~,b] = sort(ts_syn_back(:,1));
    
    ts_syn_back = ts_syn_back(b,:);
    ts_syn= [ts_syn;ts_syn_back];
    kp = abs(diff(ts_syn(:,1))-diff(ts_syn(:,3)))<1;
    ts_syn = ts_syn(kp,:);
    
    %delete points that are way off
    diff_expect  = abs(ts_syn(:,3) - ts_syn(:,1));
    kp = diff_expect < 10*median(diff_expect);
    ts_syn = ts_syn(kp,:);
    
    
    IntanTime_tmp= ts_syn(:,1);
    VideoFrameNumber_tmp = ts_syn(:,2);
    ExpectedVideoTime_tmp = ts_syn(:,3);
    
    %get the intan time for each video frame
    [VideoFrameNumber_tmp,b] = sort(VideoFrameNumber_tmp);
    
    IntanTime_tmp = IntanTime_tmp(b);
    ExpectedVideoTime_tmp = ExpectedVideoTime_tmp(b);
    
    
    %delete any repeitions
    [U,~,X] = unique(VideoFrameNumber_tmp);
    [N,E] = hist(X,1:max(X));
    C = arrayfun(@(x)find(X==x),E(N>1),'uni',0);
    kp = cell2mat(cellfun(@(a) a(2:end),C,'uni',0)');
     IntanTime_tmp(kp) =[];
    VideoFrameNumber_tmp(kp) =[];
    ExpectedVideoTime_tmp(kp) =[];
 
    IntanTime{v} = IntanTime_tmp;
    VideoFrameNumber{v} = VideoFrameNumber_tmp;
    ExpectedVideoTime{v} = ExpectedVideoTime_tmp;
    VideoNumber{v} = v*ones(length(IntanTime_tmp),1);
    %ok = fit(IntanTime, VideoFrameNumber,'poly1');
    %pt =feval(ok,IntanTime);
    
end

IntanTime = cell2mat(IntanTime');
VideoFrameNumber = cell2mat(VideoFrameNumber');
ExpectedVideoTime = cell2mat(ExpectedVideoTime');
VideoNumber = cell2mat(VideoNumber');

outTable = table(IntanTime,VideoFrameNumber,ExpectedVideoTime,VideoNumber);

outfil = [outdir filesep 'ts_syn.mat'];
save(outfil,'outTable','fils')
end