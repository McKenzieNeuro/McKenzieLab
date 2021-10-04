function [whl,in,threshF] = sm_GetLEDfromAVI(file,outname,varargin)

% This code tests an approximate median filter for separating out the red
% and blue leds mounted on top of the subjects head
% I was using the DataMax system, which was synced using a red LED. Hence,
% I am tracking the x,y position of the animal as well as when the sync
% light comes on.


plotIt = false;
readerobj  = VideoReader(file);
width      = readerobj.Width;
height     = readerobj.Height;
   out = [];
  in = [];
    threshF = [];
if length(varargin) ==1
    in = varargin{1};
    
elseif length(varargin) ==2
     in = varargin{1};
    threshF = varargin{2};
    elseif length(varargin) ==3
    
         in = varargin{1};
          threshF = varargin{2};
          
           elseif length(varargin) ==4
      
         in = varargin{1};
          threshF = varargin{2};
          plotIt = varargin{4};
          
end

% Initial frame
Fint       = 1;
x=[];y=[];
% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);

% Initialize background as a grayscale image of the first frame
bg_bw     = rgb2gray(read(readerobj,Fint));

% Initialize foreground image and difference image
fg          = zeros(size(bg_bw));
fr_diff     = zeros(size(bg_bw));

% Initialize color mask
mask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');

% Initialize whl matrix
whl = zeros(readerobj.NumberOfFrames,4);

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
in = repmat(reshape(in,[size(fr,1) size(fr,2)]),[1 1 3]);

end


%set threshold
try
    fr    = read(readerobj,round(readerobj.NumberOfFrames/2));
end
fr_bw = rgb2gray(fr);
label         = repmat(logical(fr_bw>0),[1 1 3]) & in;
mask(label)   = fr(label);
mask(~label)  = 0;


if isempty(threshF)
hh = figure;

uicontrol('Style', 'slider',...
    'Min',0,'Max',500,'Value',40,...
    'Position', [400 20 100 20],...
    'Callback', {@replot,fr_bw,hh});




%     keyboard
%%% Label Color Mask

%%% Find centroid of remaining pixels %%%
bw_mask = rgb2gray(mask);

subplot(2,1,1);imshow(uint8(fr));title('If looks wrong, Ctrl+C and change threshold !!')
ax = subplot(2,1,2);

imshow(uint8(bw_mask));title('Mask')

uiwait(hh)

end



for i = Fint:readerobj.NumberOfFrames
    %fprintf('%i',i);
    
    % Access frame of interest - if error (mostly for the last frames,
    % don't know why), will analyze previous frame...
    try
        fr    = read(readerobj,i);
    end
    % convert frame to grayscale
    fr_bw = rgb2gray(fr);
    %     keyboard
    %%% Label Color Mask
    label         = repmat(logical(fr_bw>threshF),[1 1 3]) & in ;
    mask(label)   = fr(label);
    mask(~label)  = 0;
    
    %%% Find centroid of remaining pixels %%%
    bw_mask = rgb2gray(mask);
    [L,Nl]  = bwlabel(double(bw_mask));
    ids     = [];
    
    % Initialize red and blue LEDS as missing (=> -1)
    Rr = [-1 -1];
    Br = [-1 -1];
    
    % Update centroid if connected components are found
    if Nl > 0
        a = regionprops(L,'PixelList');
        mRt = [];mBt = [];
        for ii = 1:Nl
            % Skip single pixels or large objects
            if size(a(ii).PixelList,1)<2 || size(a(ii).PixelList,1)>5000
                continue
            end
            
            % Access color information and calculate means
            R  = mask(a(ii).PixelList(:,2),a(ii).PixelList(:,1),1);
            mRt(ii) = mean(R(R>0));
            B  = mask(a(ii).PixelList(:,2),a(ii).PixelList(:,1),3);
            mBt(ii) = mean(B(B>0));
            
        end
        [mR,bR] = max(mRt);
        [mB,bB] = max(mBt);
        
        if (~isempty(mR) && ~isempty(mB)) 
            
            Rr  = round(mean(a(bR).PixelList,1));
        end
        
        
        if  (~isempty(mR) && ~isempty(mB)) 
            
            Br  = round(mean(a(bB).PixelList,1));
            
        end
        
    end
    
   
    
    whl(i,:) = [Rr(1),Rr(2),Br(1),Br(2)];
    if plotIt
    % End processing time
    if mod(i,500)==0
        h = waitbar(i/readerobj.NumberOfFrames);
        if 1
            ixr = whl(i-499:i,1);
            ixr = ixr>-1;
            ixb = whl(i-499:i,3);
            ixb = ixb>-1;
            figure(1),clf,
            subplot(3,1,1);imshow(uint8(fr));title('If looks wrong, Ctrl+C and change threshold !!')
            subplot(3,1,2),imshow(uint8(bw_mask));title('Mask')
            subplot(3,1,3)
            plot(whl(ixr,1),whl(ixr,2),'o','color','r')
            hold on
            plot(whl(ixb,3),whl(ixb,4),'o','color','b')
            fprintf('Detection of Red LED failed %i times (%i times for the blue LED) in the last 100 frames\n',sum(~ixr),sum(~ixb));
        end
    elseif mod(i,1000)==0
        h = waitbar(i/readerobj.NumberOfFrames);
    end
    end
    
    
end
save(outname,'whl');


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
