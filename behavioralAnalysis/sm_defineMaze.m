function in = sm_defineMaze(file,varargin)

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


%%


    function ImageClickCallback ( objectHandle , eventData )
        axesHandle  = get(objectHandle,'Parent');
        coordinates = get(axesHandle,'CurrentPoint');
        coordinates = coordinates(1,1:2);
        x=[x;coordinates(1)];
        y=[y;coordinates(2)];
        text(coordinates(1),coordinates(2),num2str(size(x,1)),'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
    end

    function replot(hObj,event,occMap,h)
        set(0,'currentfigure',hh)
        threshF = get(hObj,'Value');
        ax = subplot(2,1,2);
        
        label         = repmat(logical(fr_bw>threshF),[1 1 3]) & in ;
        mask(label)   = fr(label);
        mask(~label)  = 0;
        
        %%% Find centroid of remaining pixels %%%
        bw_mask = rgb2gray(mask);
        
        
        imshow(uint8(bw_mask));title('Mask');
        
        set(ax,'ButtonDownFcn',@ImageClickCallback);
    end



    function defineMaze ( objectHandle , eventData )
        close(j)
    end

end
