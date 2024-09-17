function on = sm_getStimTimes(fname,ch,varargin)

p = inputParser;



addParameter(p,'nChannels',32,@isnumeric)
addParameter(p,'Fs',20000,@isnumeric)
addParameter(p,'saveEventFile',true,@islogical)



parse(p,varargin{:})


nChannels = p.Results.nChannels;
Fs = p.Results.Fs;
saveEventFile = p.Results.saveEventFile;

dat = LoadBinary(fname,'channels',ch,'nChannels',nChannels,'frequency',Fs);
events.time = [];
events.description  =[];
for i = 1:size(dat,2)


    on{i} = find(diff(dat(:,i))>0)/Fs;
    if saveEventFile

        events.time = [events.time;on{i}];
        events.description = [events.description ;repmat({['stim on: ch ' num2str(ch(i))]},length(on{i}),1)];


    end
end
[  events.time,b] = sort(events.time);
events.description = events.description(b);
dirN = fileparts(fname);
SaveEvents([dirN filesep 'pulseTimes.evt.sti'],events)

end