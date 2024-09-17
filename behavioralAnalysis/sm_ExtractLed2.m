function [whl,in,threshF] = sm_ExtractLed2(file,varargin)

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

if ~exist([vidDir filesep 'TTL_pulse.mat'])
    [ups,dwns]  = sm_getDigitalin(vidDir,'digitalin.dat',30000,16);
else
    load([vidDir filesep 'TTL_pulse.mat'])
end


fs = readerobj.FrameRate;
threshF = prctile(diff(whl),97);
up_vid = find(diff([0;[0;diff(whl)>threshF]])>0);
dwn_vid = find(diff([0;[0;diff(whl)<-threshF]])>0);

dt = median(dwns{1}-ups{1});
dt_fr = median(round(dt*fs));

kp = abs((up_vid+dt_fr)-bestmatch(up_vid+dt_fr,dwn_vid))<5;
up_vid = up_vid(kp);



max_dt  =max(diff(ups{1}))+.2;
%%
% now find matches
minJ = 1;ts_syn =[];
it=1;
idx=1;

allfound = false;

while ~allfound
    
    if it<length(ups{1})-3
        dt3 = diff(ups{1}(it:(it+3)));
        dt3_vid = diff(up_vid(minJ:minJ+3))/fs;
        
        if all(abs(dt3 - dt3_vid)<.15)
            
            ts_syn = [ts_syn;ups{1}(it) up_vid(minJ) up_vid(minJ)/fs minJ];
            
            
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
                    for ii = it:length(ups{1})-4
                        dt3 = diff(ups{1}(ii:(ii+4)));
                        if all(abs(dt3 - dt3_vid)<.15)
                            
                            % check if it is on regression line
                            
                            if size(ts_syn,1)>=3
                            kp=~isnan(ts_syn(:,2));
                            ok = fit(ts_syn(kp,1) , ts_syn(kp,2),'poly1');
                            pt =feval(ok,ups{1}(ii));
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
                    
                    if ii==length(ups{1})-4 && jj ==length(up_vid)-4
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
if size(ups{1},1) ~= size(ts_syn,1)
    
    
    
    % now find matches
    minJ = length(up_vid)-3;
    it=length(ups{1})-3;
    idx=1;
    ts_syn_back =[];
    allfound = false;
    
    
    while ~allfound
        
        if it>0
            dt3 = diff(ups{1}(it:(it+3)));
            dt3_vid = diff(up_vid(minJ:minJ+3))/fs;
            
            if all(abs(dt3 - dt3_vid)<.15)
                
                ts_syn_back = [ts_syn_back;ups{1}(it) up_vid(minJ) up_vid(minJ)/fs minJ];
                
                
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
                            dt3 = diff(ups{1}(ii:(ii+3)));
                            if all(abs(dt3 - dt3_vid)<.15)
                                
                                % check if it is on regression line
                                
                                if size(ts_syn_back,1)>=3
                                    kp=~isnan(ts_syn_back(:,2));
                                    ok = fit(ts_syn_back(kp,1) , ts_syn_back(kp,2),'poly1');
                                    pt =feval(ok,ups{1}(ii));
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
outfil = [vidDir filesep 'ts_syn.mat'];
save(outfil,'ts_syn','in')
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
