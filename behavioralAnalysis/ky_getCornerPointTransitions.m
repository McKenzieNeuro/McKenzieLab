function ky_getCornerPointTransitions
% load java excel libraries for linux system
% javaaddpath('poi_library/poi-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('poi_library/dom4j-1.6.1.jar');
% javaaddpath('poi_library/stax-api-1.0.1.jar');



prompt = {'Enter Event IDs (comma separated)'};
dlg_title = 'Event ID';
num_lines = 1;
def = {''};
ID = [];
data = [];
tbl = [];
c = [];
f = [];
t = [];
IDs = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(IDs)
C = strsplit(IDs{1},',');
C = C(:);

%%

[fname,dirname]=uigetfile('*.mp4');
fname=[dirname fname];

if isstr(fname)
    if(exist('fname','var')~=1); vid = [];
    elseif(ischar(fname) && exist(fname,'file')==2)
        vid = VideoReader(fname);
    elseif(strcmp(class(fname),'VideoReader'))
        vid = fname;
    elseif(isa(fname,'uint8') && size(fname,3)==3)
        vid = []; staticframe = fname;
    else
        error('retrackLAB:InvalidArgument',...
            '''fname'' input argument should be filename or ''VideoReader'' class.');
    end
    
    
    
    
    % Process frame or timestamp specified as input:
    if(exist('ts','var')==1 && ~isempty(ts))
        frame = find(AVIts>=ts,1,'first');
    end
    if(exist('frame','var')~=1 || isempty(frame))
        frame = 1;
    end
    
    
    % Initialize frame rate.
    rate = 1;
    skip = 1;
    dt = round(1/(vid.FrameRate*rate)*1000)/1000;
    presses=[];
    
    % Context geometry state
    context = struct([]);
    contextCounter = 0;

    % Create timer to drive the video
    t = timer('TimerFcn',@timerfcn,'Period',dt,'ExecutionMode','fixedDelay',...
        'UserData',frame,'StartFcn',@startfcn,'StopFcn',@stopfcn);
    
    
    %%
    c = figure('CloseRequestFcn',@(~,~) closeButton([],[]));
    
    scr =get(0,'MonitorPositions');
         scr = scr(1,:);
    set(c,'position', [20   scr(4)-100*length(C)   215   70*length(C)],'name','No ID','NumberTitle','off')
    for i = 1:length(C)
        IDb = uicontrol('Parent',c,'Style','pushbutton','String',C{i},...
            'Units','normalized','Position',[0 1-i/length(C) 1 1/(length(C)+1) ],'Units','normalized',...
            'Tag','ID','Callback',@setID);
        
    end
    
    %%
    % Create the video playback GUI
    f = figure('KeyPressFcn',@(x,y) clockFrame(x,y),'CloseRequestFcn',@closeButton);
    
    
    p = get(f,'Position');
    set(f,'Position',[p(1) p(2) 860 600]);
    movegui(f,'center');
    
    a = axes('Parent',f,'Units','Pixels','Position',[10 110 640 480],...
        'DataAspectRatio',[1 1 1],'XLim',[0.5 640.5],'YLim',[.5 480.5],...
        'Box','on','XTick',[],'YTick',[],'LooseInset',[0 0 0 0],...
        'Units','normalized');
    if(~isempty(vid))
        im = image('XData',[1 640],'YData',[480 1],'CData',read(vid,frame),'Parent',a);
    elseif(exist('staticframe','var')==1)
        im = image('XData',[1 640],'YData',[480 1],'CData',staticframe,'Parent',a);
    end
    
    ui_frame = uicontrol('Parent',f,'Style','edit','String',num2str(frame),...
        'Units','pixels','Position',[ 15 70 100 30],'Units','normalized',...
        'Tag','frame','Callback',@changeval);
    
    ui_rate = uicontrol('Parent',f,'Style','edit','String',num2str(rate),...
        'Units','pixels','Position',[399 70 100 30],'Units','normalized',...
        'Tag','rate','Callback',@changeval);
    ui_save = uicontrol('Parent',f,'Style','pushbutton','String','Save Clicks',...
        'Units','pixels','Position',[527 70 100 30],'Units','normalized',...
        'Callback',@saveAs);
    
    uicontrol('Parent',f,'Style','pushbutton','String','Prev',...
        'Units','pixels','Position',[ 15 10 100 50],'Units','normalized',...
        'Callback',@prev);
    ui_play = uicontrol('Parent',f,'Style','pushbutton','String','Play',...
        'Units','pixels','Position',[137 10 100 50],'Units','normalized',...
        'Callback',@play);
    uicontrol('Parent',f,'Style','pushbutton','String','Next',...
        'Units','pixels','Position',[259 10 100 50],'Units','normalized',...
        'Callback',@next);
    uicontrol('Parent',f,'Style','pushbutton','String','Slower',...
        'Units','pixels','Position',[381 10 100 50],'Units','normalized',...
        'Callback',@slow);
    uicontrol('Parent',f,'Style','pushbutton','String','Faster',...
        'Units','pixels','Position',[503 10 100 50],'Units','normalized',...
        'Callback',@fast);
    
    uicontrol('Parent',f,'Style','pushbutton','String','Undo',...
        'Units','pixels','Position',[625 10 100 50],'Units','normalized',...
        'Callback',@undo);
    
    uicontrol('Parent',f,'Style','pushbutton','String','Show data',...
        'Units','pixels','Position',[747 10 100 50],'Units','normalized',...
        'Callback',@showData);

    uicontrol('Parent',f,'Style','pushbutton','String','Mark perimeter',...
        'Units','pixels','Position',[747 70 100 50],'Units','normalized',...
        'Callback',@addNewEvent);
    
    clicksClocked = text(150,-40,'Total Clicks:');
    LastClick = text(270,-40,'Last Click:');
    
    
end
end

    function dispframe(n,c)
        
        if(n < vid.NumberOfFrames && n > 0)
            if(~isempty(vid))
                try
                    set(im,'CData',read(vid,n));
                catch
                end
            end
            
            set(ui_frame,'String',num2str(n));
            
            t.UserData = n;
            
            
        end
    end

    function startfcn(~,~)
        set(ui_play,'Enable','on','String','Pause');
    end

    function stopfcn(~,~)
        set(ui_play,'Enable','on','String','Play');
    end

    function timerfcn(obj,~)
        frame = obj.UserData + skip;
        if(frame > vid.NumberOfFrames || frame <= 0); stop(obj); end
        dispframe(frame);
        drawnow;
    end
    function setID(obj,test)
        ID =  get(obj,'String');
        ID = find(cellfun(@(a) strcmp(a,ID),C));
        set(get(obj,'Parent'),'name',C{ID},'NumberTitle','off')
        figure(f)
    end
    function play(~,~)
        if(strcmp(t.Running,'on'))
            set(ui_play,'Enable','off','String','Play');
            stop(t);
        else
            set(ui_play,'Enable','off','String','Pause');
            start(t);
        end
        ReleaseFocus(f)
    end

    function saveAs(source,~)
        if(strcmp(t.Running,'on'))
            set(ui_play,'Enable','off','String','Play');
            stop(t);
        end
        closeButton(source,1);
    end



    function prev(~,~)
        stop(t);
        frame = t.UserData - 1;
        dispframe(frame);
        ReleaseFocus(f)
    end

    function next(~,~)
        stop(t);
        frame = t.UserData + 1;
        dispframe(frame);
        ReleaseFocus(f)
    end

    function changerate(r)
        if(r == 0); r = rate; end
        if(abs(r)>1);
            rate = round(r);
            %rate = 1;
            dt = round(1/vid.FrameRate*1000)/1000;
            skip = rate;
        else
            rate = r;
            dt = round(1/(vid.FrameRate*abs(rate))*1000)/1000;
            skip = sign(rate);
        end
        
        if(strcmp(t.Running,'on'))
            stop(t);
            set(t,'Period',dt);
            start(t);
        else set(t,'Period',dt);
        end
        set(ui_rate,'String',num2str(rate));
    end

    function slow(~,~)
        changerate(sign(rate)*2^(round(log2(abs(rate)))-1))
        ReleaseFocus(f)
    end

    function fast(~,~)
        changerate(sign(rate)*2^(round(log2(abs(rate)))+1));
        ReleaseFocus(f)
    end

    function undo(~,~)
        if ~isempty(presses)
            presses=presses(1:end-1,:);
            
            if ~isempty((presses))
                set(clicksClocked,'string',['Total Clicks:' num2str(size(presses,1))]);
                set(LastClick,'string',['Last click:' num2str(presses(end,1))]);
            else
                set(clicksClocked,'string',['Total Clicks: 0']);
                set(LastClick,'string',['Last click: ']);
                
            end
        end
        
        
    end

    function clockFrame(x,y)
        
        if isempty(ID)
            warndlg('SELECT YOUR KEY')
            
            
            if(strcmp(t.Running,'on'))
                set(ui_play,'Enable','off','String','Play');
                stop(t);
            end
            
            
        else
            
            if isempty(presses) || ~ismember(frame,presses(:,1))
                
                
                presses=[presses;frame ID];
                
                
                set(clicksClocked,'string',['Total Clicks:' num2str(size(presses,1))]);
                set(LastClick,'string',['Last click:' num2str(presses(end,1))]);
                
                
                if(strcmp(t.Running,'on'))
                    set(ui_play,'Enable','off','String','Play');
                    stop(t);
                end
            elseif strcmp(t.Running,'off')
                set(ui_play,'Enable','off','String','Pause');
                start(t);
            end
            
           
            ReleaseFocus(f)
        end
    end


    function changeval(obj,~)
        val = str2double(get(obj,'String'));
        switch get(obj,'Tag')
            case 'frame'
                dispframe(round(val));
            case 'coord'
                n = find(COMmap==val,1,'first');
                dispframe(n,val);
            case 'ts'
                val = find(AVIts>=val,1,'first');
                dispframe(val);
            case 'rate'
                changerate(val);
        end
    end

    function ReleaseFocus(fig)
        if ~exist('fig','var') || ~isvalid(fig) % Guard callback that uses figure
            return
        end
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'off');
        drawnow;
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'on');
    end


    function showData(obj,~)
        fig = figure('CloseRequestFcn',@closeTable);
        tbl = axes('units', 'normalized','position', [0 0 .8 .8]);
        columninfo.titles={'ID','ts'};
        columninfo.formats = {'%4.6g','%4.6g','%4.6g','%4.6g', '%4.6g'};
        columninfo.weight =      [ 1, 1];
        columninfo.multipliers = [ 1, 1];
        columninfo.isEditable =  [ 1, 1];
        columninfo.isNumeric =   [ 0  1];
        columninfo.withCheck = false; % optional to put checkboxes along left side
        columninfo.chkLabel = 'Use'; % optional col header for checkboxes
        rowHeight = 16;
        gFont.size=9;
        gFont.name='Helvetica';
        if any(presses)
            figData = mltable(fig, tbl, 'CreateTable', columninfo, rowHeight, [C(presses(:,2)) num2cell(presses(:,1),2)], gFont);
            uiwait(fig)
        end
        ReleaseFocus(f)
        
        function closeTable(obj,~)
            ok = get(tbl,'userdata');
            if ~isempty(ok)
                ok = ok.data;
                kp = ismember(presses(:,1),cell2mat(ok(:,2)));
                presses = presses(kp,:);
                
                if ~isempty((presses))
                    set(clicksClocked,'string',['Total Clicks:' num2str(size(presses,1))]);
                    set(LastClick,'string',['Last click:' num2str(presses(end,1))]);
                else
                    set(clicksClocked,'string',['Total Clicks: 0']);
                    set(LastClick,'string',['Last click: ']);
                    
                end
                
                
            end
            
           
            delete(obj)
            
        end
    end

    function addNewEvent(~,~)
        % Stop playback
        if exist('t','var') && isvalid(t) && strcmp(t.Running,'on')
            stop(t);
        end

        contextCounter = contextCounter + 1;
        k = contextCounter;

        % Ask user for real-world dimensions
        prompt = {'Length:','Width:','Units (cm or in):','Bin size (cm):'};
        def = {'40', '40', 'cm', '1'};
        answ = inputdlg(prompt,'Context Geometry',1,def);
        if isempty(answ)
            contextCounter = contextCounter - 1;
            return;
        end

        L = str2double(answ{1});
        W = str2double(answ{2});
        units = lower(answ{3});
        bin_cm = str2double(answ{4});

        if strcmp(units,'in')
            L_cm = L * 2.54;
            W_cm = W * 2.54;
        else
            L_cm = L;
            W_cm = W;
        end

        % Grab current frame
        frameNum = t.UserData;
        frameImg = read(vid, frameNum);

        % Corner selection UI
        figCP = figure('Name','Select Context Corners','NumberTitle','off');
        imshow(frameImg);hold on;

        h = size(frameImg,1);

        % Visual axis lock
        quiver(20,h-20,80,0,'r','LineWidth',2);
        text(110,h-25,'Length (X)','Color','r','FontWeight','bold');

        quiver(20,h-20,0,-80,'b','LineWidth',2);
        text(25,h-110,'Width (Y)','Color','b','FontWeight','bold');

        title({ ...
            'Click 4 corners CLOCKWISE' ...
            'Start at BOTTOM LEFT' ...
            'Length → X | Width ↓ Y' ...
            });

        [x,y] = getpts;
        if isvalid(figCP); close(figCP); end

        if numel(x) ~=4
            warndlg('You must select exactly 4 corners.');
            contextCounter = contextCounter - 1;
            return;
        end

        imagePts = enforceClockwise([x y]);

        % ---- World-space corners ----
        worldPts = [
            0     0
            L_cm  0
            L_cm  W_cm
            0     W_cm
            ];

        % ---- Fit geometric transform ----
        tform = fitgeotrans(imagePts, worldPts, 'projective');

        % ---- Pixel/cm ratio ----
        pixL = norm(imagePts(2,:) - imagePts(1,:));
        pixW = norm(imagePts(3,:) - imagePts(2,:));
        px_per_cm = mean([pixL/L_cm, pixW/W_cm]);

        % ---- Warn on geometry jump ----
        if k > 1
            prev = context(k-1).px_per_cm;
            jump = abs(px_per_cm - prev) / prev;
            if jump > 0.15
                warndlg(sprintf('Pixel/cm changed by %.1f%%',jump*100), ...
                    'Geometry Warning');
            end
        end

        % ---- Warp preview ----
        ref = imref2d([ceil(W_cm/bin_cm), ceil(L_cm/bin_cm)], ...
            [0 L_cm], [0 W_cm]);
        warped = imwarp(frameImg, tform, 'OutputView', ref);

        figPrev = figure('Name','Context Preview','NumberTitle','off');
        imshow(warped); grid on;
        set(gca,'XTick',0:bin_cm:L_cm,'YTick',0:bin_cm:W_cm)
        xlabel('Length (cm)');
        ylabel('Width (cm)');
        title(sprintf('Preview: %.1f × %.1f cm',L_cm,W_cm));

        % ---- Mask + meshgrid ----
        [X,Y] = meshgrid(1:size(frameImg,2),1:size(frameImg,1));
        mask = poly2mask(imagePts(:,1),imagePts(:,2), ...
            size(frameImg,1),size(frameImg,2));
        pixelCoords = [X(mask) Y(mask)];

        % ---- Store context ----
        context(k).eventID      = ID;
        context(k).startFrame  = frameNum;
        context(k).length_cm   = L_cm;
        context(k).width_cm    = W_cm;
        context(k).bin_cm      = bin_cm;

        context(k).imagePts    = imagePts;
        context(k).worldPts    = worldPts;
        context(k).tform       = tform;
        context(k).px_per_cm   = px_per_cm;

        context(k).mask        = mask;
        context(k).pixelCoords = pixelCoords;
        context(k).name        = sprintf('Context_%gx%gcm',L_cm,W_cm);

        disp(['Added ' context(k).name ' at frame ' num2str(frameNum)]);
    end

    function closeButton(~,saved)
        % Safe GUI shutdown (can be called multiple times)

        % ---- Stop & delete timer safely ----
        if exist('t','var') && isa(t,'timer') && isvalid(t)
            try
                if strcmp(t.Running,'on')
                    stop(t);
                end
                delete(t);
            catch
            end
        end

        % ---- Save logic ----
        if ~isempty(presses) || ~isempty(context)
            if nargin < 2 || isempty(saved)
                choice = questdlg('Would you like to save?', ...
                    'Save Data', 'Yes','No','Yes');
            else
                choice = 'Yes';
            end

            if strcmp(choice,'Yes')
                [FileName,dname] = uiputfile('*.mat');
                if ischar(FileName)
                    fout = fullfile(dname,FileName);
                    save(fout,'presses','context')
                    
                end
            end
        end

        % ---- Close figures safely ----
        if exist('c','var') && isvalid(c)
            delete(c)
        end
        if exist('f','var') && isvalid(f)
            delete(f)
        end
    end

    function pts = enforceClockwise(pts)
        ctr = mean(pts,1);
        ang = atan2(pts(:,2)-ctr(2), pts(:,1)-ctr(1));
        [~,idx] = sort(ang);
        pts = pts(idx,:);
    end




end
