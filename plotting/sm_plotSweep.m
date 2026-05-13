function sm_plotSweep(fname,time,duration,varargin)


set(0, 'DefaultFigureRenderer', 'painters');


p = inputParser;
addParameter(p,'channels',[],@isnumeric);
addParameter(p,'shank',[],@isnumeric);
addParameter(p,'scaleFactor',.195,@isnumeric);
addParameter(p,'fout',[],@isstring);
addParameter(p,'color',[],@isnumeric);

parse(p,varargin{:});

channels = p.Results.channels;
shank = p.Results.shank;
scaleFactor = p.Results.scaleFactor;
fout = p.Results.fout;
color = p.Results.color;

gap  = 700;
if ~exist(fname)
    error('file does not exist')
end

fxml = strrep(fname,'.dat','.xml');

if exist(fxml)
    
    xml = LoadXml(fxml);
else
    
    error('make xml')
    
end

if isempty(channels) & isempty(shank)
    channels = {0:xml.nChannels-1};
elseif ~isempty(shank)
    
    if ~isempty(channels)
        warning('shank assignment overrides channel')
    end
    if max(shank) > length(xml.AnatGrps)
        error('too many shanks')
        
    end
    
    channels = {xml.AnatGrps(shank).Channels};
    
elseif ~isempty(channels)
    
    channels = {channels};
    
end

if ~isempty(color)
    if size(color,1)>1
        error('only 1 color supported')
    end
    color = {color};
else
    
    color = linspecer(length(channels),'jet');
end



if isempty(fout)
    fout = strrep(fname,'.dat',['_sweep' num2str(round(time)) '.ps']);
    
    
end

% load data
h = figure;

idx = 1;
for sh = 1:length(channels)
   d = LoadBinary(fname,'nchannels',xml.nChannels,'channels',channels{sh}+1,'frequency', ...
       xml.SampleRate,'start',time,'duration',duration);
    d = double(d) * scaleFactor;
    
    ts = (1:size(d,1))/xml.SampleRate + time;
    
    for i = 1:size(d,2)
        plot(ts,d(:,i) - (idx-1)*gap,'color',color{sh},'linewidth',2)
        hold on
        idx = idx+1;
    end
end


plot(7*[range(ts)/8 range(ts)/8]+time,-[(idx-1)*gap+2000 (idx-1)*gap+4000],'k','linewidth',1)
text(7*range(ts)/7.9 +time,-[(idx-1)*gap+3000],'2mV')
set(gca,'ytick',[],'fontsize',16)
ylim([-((idx-1)*gap+6000) 2000])
%set(gcf,'position',[  1          41        1920        1083])
set(h, 'PaperUnits', 'inches');
set(h, 'PaperOrientation', 'landscape');   % switch to landscape mode
set(h, 'PaperPositionMode', 'manual');

set(h, 'PaperSize', [11 8.5]);             % width x height for US Letter (landscape)
set(h, 'PaperPosition', [0 0 11 8.5]);
print(h, '-dpsc2',fout );

filenamepdf = strrep(fout,'ps','pdf');

 sm_ps2pdf(fout,filenamepdf,[])

