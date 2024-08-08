function sm_getArenaEdges(dirN)
% load java excel libraries for linux system
% javaaddpath('poi_library/poi-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('poi_library/dom4j-1.6.1.jar');
% javaaddpath('poi_library/stax-api-1.0.1.jar');

% cd('R:\McKenzieLab\DANEHippocampalResponse\DA2h2\Novel Environment\DA2h2-210817-141645')

% cd('R:\McKenzieLab\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220720-123755')

% cd('R:\McKenzieLab\ASommer\FP experiments DA-NE\Novel Environment\Full Experiment\Day 1\DA2h9-230130-095543')
cd(dirN)
prompt = {'Enter Context IDs (comma separated)'};
dlg_title = 'Context ID';
num_lines = 1;
def = {''};
ID = [];

data = [];
tbl = [];
IDs = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(IDs)
    C = strsplit(IDs{1},',');
    
    
    %%
    
    [fname,dirname]=uigetfile('*.avi');
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
        % Create timer to drive the video
        t = timer('TimerFcn',@timerfcn,'Period',dt,'ExecutionMode','fixedDelay',...
            'UserData',frame,'StartFcn',@startfcn,'StopFcn',@stopfcn);
        
        
        %%
        c = figure;
        
        scr =get(0,'MonitorPositions');
        set(c,'position', [20   scr(4)/2+500  215   70*length(C)],'name','No ID','NumberTitle','off')
        for i = 1:length(C)
            IDb = uicontrol('Parent',c,'Style','pushbutton','String',C{i},...
                'Units','normalized','Position',[0 1-i/length(C) 1 1/(length(C)+1) ],'Units','normalized',...
                'Tag','ID','Callback',@setID);
            
        end
        coords = cell(length(c),1);
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
        uicontrol('Parent',f,'Style','pushbutton','String','Next',...
        'Units','pixels','Position',[259 10 100 50],'Units','normalized',...
        'Callback',@next);
        uicontrol('Parent',f,'Style','pushbutton','String','Prev',...
            'Units','pixels','Position',[ 15 10 100 50],'Units','normalized',...
            'Callback',@prev);
        ui_play = uicontrol('Parent',f,'Style','pushbutton','String','Play',...
            'Units','pixels','Position',[137 10 100 50],'Units','normalized',...
            'Callback',@play);
        
        uicontrol('Parent',f,'Style','pushbutton','String','Slower',...
            'Units','pixels','Position',[381 10 100 50],'Units','normalized',...
            'Callback',@slow);
        uicontrol('Parent',f,'Style','pushbutton','String','Faster',...
            'Units','pixels','Position',[503 10 100 50],'Units','normalized',...
            'Callback',@fast);
        
        uicontrol('Parent',f,'Style','pushbutton','String','LabelCoords',...
            'Units','pixels','Position',[625 10 100 50],'Units','normalized',...
            'Callback',@getCoords);
        
        
        
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

 function getCoords(obj,test)
       coords{ID} = ginput(4); 
 end


    function play(~,~)
        if(strcmp(t.Running,'on'));
            set(ui_play,'Enable','off','String','Play');
            stop(t);
        else
            set(ui_play,'Enable','off','String','Pause');
            start(t);
        end
        ReleaseFocus(f)
    end

    function saveAs(source,~)
        if(strcmp(t.Running,'on'));
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
        if(abs(r)>1)
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
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'off');
        drawnow;
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'on');
    end

    function closeButton(source,saved) %#ok<VANUS>
        %runs when GUI is shut or when user exits, prompts for saving data
        %temporary file only deleted if user explicitly says not to save
        
        %stop timers
        
        if isvalid(t) && strcmp(t.running,'on')
            stop(t)
        end
        
        if ~ all(cellfun(@isempty,coords)==1)
            %prompt user to save
            if isempty(saved)
                choice = questdlg('Would you like to save?', ...
                    'Save Data', 'Yes','No','Yes');
            else
                choice='Yes';
            end
            switch choice
                case 'Yes'
                    
                    
                    [FileName,dname] = uiputfile('*.mat','');
                    
                    
                    
                    if ischar(FileName)
                        
                        contextEdges.contextName = C;
                        contextEdges.edges = coords;
                        save([dname filesep FileName],'contextEdges');
                        
                        
                        
                        
                        
                        
                    else
                        warndlg('No coords clocked, no file saved','Enjoy yourself!')
                    end
                    
                    
                    
                    
            end
        else
            warndlg('No clicks clocked, no file saved','Enjoy yourself!')
        end
        close(c)
        delete(f)
    end



end
