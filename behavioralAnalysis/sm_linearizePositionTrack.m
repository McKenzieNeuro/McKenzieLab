function  [len,pos,xp,yp,del] =sm_linearizePositionTrack(X,Y)
close all
X(X==-1) = nan;
Y(Y==-1) =nan;
del = X==0 & Y==0;
len=[];
pos=[];

x=[];
y=[];
%visualize the maze without a real picture
binX = 0:max(X);
binY = 0:max(Y);
occMap = histcn([X Y],binX,binY);
occMap = SmoothMat(occMap, [5 5], 3);
h = figure;

occMap(occMap==0) =nan;
%define 20 points around the maze
ax = imagesc(occMap'>.07);
uicontrol('Style', 'slider',...
    'Min',0,'Max',nanmedian(occMap(:))*20,'Value',.07,...
    'Position', [400 20 100 20],...
    'Callback', {@replot,occMap,h});

set(ax,'ButtonDownFcn',@ImageClickCallback);

uicontrol('Style', 'pushbutton', 'String', 'defineMaze',...
    'Position', [20 20 100 20],...
    'Callback', @defineMaze);

waitfor(h)


    function  defineMaze( objectHandle , eventData )
        
        tolerance = 500;
        %find closest point on the maze
        [qx,qy] = getPoints([x;x(1)]',[y;y(1)]',2600);
        trackpos = [qx' qy'];
        [trackcoordInd,mindist]=dsearchn(trackpos,[X Y]);
        del = del | mindist > tolerance;
        xp = trackpos(trackcoordInd,1);
        yp = trackpos(trackcoordInd,2);
        
        
        
        %find leftmost point
        xl = min(xp); yl = min(yp);
        
        %rightmost
        xr = max(xp); yr = max(yp);
        
        cx = xp-x(1);
        cy = yp - y(1);
        
        len = hypot(cx,cy);
        len = inpaint_nans(len);
        j = figure;
           uicontrol('Style', 'pushbutton', 'String', 'redo',...
            'Position', [20 20 100 20]     );
      hold on
        ax = plot(X,Y,'x');
        
        
     
        
        hold on
        
         
      [x,y] = getpts;
        
        
       % waitfor j
        
        IN = inpolygon(X,Y,x,y);
        del = del | ~IN;
        %get rid of tracking errors
        len(del) =nan;
        
        %bin the space now in terms of distance from a fixed point
        binsize = 10;
        bins = 0:binsize:max(len)+binsize;
        [occ,pos] = histc(len,bins);
     
        
        %
    end

    function ImageClickCallback ( objectHandle , eventData )
        axesHandle  = get(objectHandle,'Parent');
        coordinates = get(axesHandle,'CurrentPoint');
        coordinates = coordinates(1,1:2);
        x=[x;coordinates(1)];
        y=[y;coordinates(2)];
        text(coordinates(1),coordinates(2),num2str(size(x,1)),'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
    end


    function drawbox ( objectHandle , eventData )
        
        clf(objectHandle,'reset')
        plot(X(1:100:end),Y(1:100:end),'x')
        hold on
        x=[];
        y=[];
      
      [x,y] = ginput(4);
       
               plot(x,y,'--mo')
      
        plot([x;x(1)],[y;y(1)],'--mo')
        
        
       
        
    end


    function replot(hObj,event,occMap,h)
        set(0,'currentfigure',h)
        thres = get(hObj,'Value');
        ax = imagesc(occMap'>thres);
        text(x,y,num2cell(1:size(x,1),2),'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
        set(ax,'ButtonDownFcn',@ImageClickCallback);
    end
end



