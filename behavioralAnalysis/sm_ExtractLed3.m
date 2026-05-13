function [whl,in,threshF,fs] = sm_ExtractLed3(file,varargin)

p = inputParser;
addParameter(p,'threshF',1,@isnumeric);
addParameter(p,'in',[],@islogical);
addParameter(p,'plotIt',true,@islogical);
addParameter(p,'inputch',1,@isnumeric); % red = ch1


parse(p,varargin{:});

threshF = p.Results.threshF;
in = p.Results.in;
plotIt = p.Results.plotIt;
inputch = p.Results.inputch;
% ts_video = interp1(outTable.VideoFrameNumber,outTable.IntanTime,1:maxFrame,'linear','extrap');

%%


fout = [file(1:end-4) '_LED.mat'];
readerobj  = VideoReader(file);
width      = readerobj.Width;
height     = readerobj.Height;
vidDir = fileparts(file);
% Initial frame
Fint       = 1;
x=[];y=[];
% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);


% Initialize color mask
mask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');

% Initialize whl matrix
whl = zeros(readerobj.NumberOfFrames,1);
fs = readerobj.FrameRate;

%get maze area
%
j=  figure;
try
    fr    = read(readerobj,1000);
end


if isempty(in)
    ax = imshow(uint8(fr));
    set(ax,'ButtonDownFcn',@ImageClickCallback);
    
    uicontrol('Style', 'pushbutton', 'String', 'defineMaze',...
        'Position', [20 20 100 20],...
        'Callback', @defineMaze);
    
    
    
    
    
    
    
    waitfor(j)
    
    
    
    
    
    [XX,YY] = meshgrid(1:size(fr,2),1:size(fr,1));
    in = inpolygon(X(:),Y(:),x,y);
    
    
end


try
    fr    = read(readerobj,round(readerobj.NumberOfFrames/2));
end





for i = Fint:readerobj.NumberOfFrames
    %fprintf('%i',i);
    
    % Access frame of interest - if error (mostly for the last frames,
    % don't know why), will analyze previous frame...
    try
        fr    = read(readerobj,i);
    end
    % convert frame to grayscale
    fr = fr(:,:,inputch);
    
    whl(i,1) = nanmean(fr(in));
    
end

save(fout,'whl','in','fs')

%%


function ImageClickCallback ( objectHandle , eventData )
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint');
coordinates = coordinates(1,1:2);
x=[x;coordinates(1)];
y=[y;coordinates(2)];
text(coordinates(1),coordinates(2),num2str(size(x,1)),'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
end



function defineMaze ( objectHandle , eventData )
close(j)
end

end
