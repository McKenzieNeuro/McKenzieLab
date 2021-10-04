function   sm_MakeEventTime(fname,evs,varargin)
% makes an event file with the start/stop/center times for either square pulse events
% or half sine events delivered to an analogin in channel.
dirName = fileparts(fname);
thres1 = 1000;
thres2 = 1000;
ch = 1:4;
events1.time =[];
events1.description =[];thres =[];
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'ch'
            ch = varargin{i+1};
            ch = ch(:)';
        case 'thres1'
            thres1 = varargin{i+1};
            case 'thres2'
            thres2 = varargin{i+1}; 
          
            
    end
    
end



pulse = [];

% get pulse/sine epoch labels from user define *.sti file
pulse_on =  unique(evs.time(cellfun(@any,regexpi(evs.description,'pulse_on'))));
pulse_off =  unique(evs.time(cellfun(@any,regexpi(evs.description,'pulse_off'))));

sine_on =  unique(evs.time(cellfun(@any,regexpi(evs.description,'sine_on'))));
sine_off =  unique(evs.time(cellfun(@any,regexpi(evs.description,'sine_off'))));


%check to see if there are matched onset and offset
if length(pulse_on) ~= length(pulse_off)
    disp('mismatch pulse onset/offset')
    
    return;
end

if length(sine_on) ~= length(sine_off)
    disp('mismatch sine onset/offset')
    
    return;
end


%loop through analogin channel (n = 4) and find pulses/sines
idx1 = 1;idx = 1;
events.time = [];
events.description = [];
%filename = getAllExtFiles(dirName,'xml',0);
%filename = filename(cellfun(@any,regexp(filename,'analogin')));
filename = [fname(1:end-3) 'xml'];
[xml, rxml] = LoadXml(filename);
for i = ch
    
    
    
    %look for either pulses or sines
    for kk = 1:2
        onset = [];
        offset = [];
        center =[];
        pulse = [];
        switch kk
            case 1
                if any(pulse_on)
                    
                    %detect pulse events
                    [onset,offset,pulse,center] = sm_getPulseTime(fname,i,'pulse','events',[pulse_on pulse_off],'thres',thres1);
                    label = 'pulse';
                end
                
            case 2
                if any(sine_on)
                    label = 'sine';
                    %detect sine events
                    [onset,offset,pulse,center] = sm_getPulseTime(fname,i,'sine','events',[sine_on sine_off],'thres',thres1,'thres1',thres2);
                end
        end
        
        
        
        if    ~isempty(center)
            
            %store all events in main struct to be saved
            onset = cellfun(@(a) a(:) ,onset,'uni',0);
            offset = cellfun(@(a) a(:) ,offset,'uni',0);
            center = cellfun(@(a) a(:) ,center,'uni',0);
            pulse = cellfun(@(a) a(:) ,pulse,'uni',0);
            onset = cell2mat(onset');
            offset = cell2mat(offset');
            center = cell2mat(center');
            pulse = cell2mat(pulse');
            for j = 1:length(onset)
                events.time(idx) =  onset(j);
                events.description{idx} =   [label '_on ' num2str(i)];
                idx = idx+1;
            end
            
            
            
            
            
            for j = 1:length(offset)
                events.time(idx) =  offset(j);
                events.description{idx} =   [label '_off ' num2str(i)];
                idx = idx+1;
            end
            
            for j = 1:length(center)
                events1.time(idx1) =  center(j);
                events1.description{idx1} =   [label '_center ' num2str(i) ' ' num2str(pulse(j))];
                idx1 = idx1+1;
            end
            
            
            
            
        end
    end
    
end


%sort events by time
[events.time,b] = sort(events.time);
%events.time = events.time+6072.185;
events.description = events.description(b);

[~,f]  =fileparts(fname);
outfile = [f '.evt.ait'];
sm_SaveEvents(outfile,events)

%sort events by time
[events1.time,b] = sort(events1.time);
%events.time = events.time+6072.185;
events1.description = events1.description(b);



outfile = [f '.evt.aip'];
sm_SaveEvents(outfile,events1)

end

