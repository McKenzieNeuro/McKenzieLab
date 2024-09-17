function step11_IIICGUI_mPatients
    close all;clc;warning('off','all')
   
    dataDir='./Data/CP_centers/';
    
    user_ID=inputdlg('Please enter your initials','Hi');
    rater=user_ID{1};
    targetDir=[pwd,'/Output/',rater,'/'];
    
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    set(f,'HitTest','off')
    set(f,'WindowButtonDownFcn',@clicks_Callback);
    set(f,'ButtonDownFcn',@clicks_Callback);
    set(f,'KeyPressFcn',@keys_Callback);
    set(f,'MenuBar','none','ToolBar','none');

    f1=figure('units','normalized','outerposition',[0 0 .5 1]);
    set(f1,'HitTest','off')
    set(f1,'KeyPressFcn',@keys_Callback);
    set(f1,'Visible','off');
    set(f1,'MenuBar','none','ToolBar','none');
    
    patterns={'Seizure','LPD','GPD','LRDA','GRDA','Other'};
    nP=length(patterns);
    shortcuts={'1','2','3','4','5','6'};

    set(0,'CurrentFigure',f)
    Ax_spec={subplot('position',[.03 .445 .26 .11]);
             subplot('position',[.03 .330 .26 .11]);
             subplot('position',[.03 .215 .26 .11]);
             subplot('position',[.03 .100 .26 .11]);
             subplot('position',[.03 .560 .26 .01])};
    Ax_EEG=subplot('position',[.32 .1 .6 .85]);
    Ax_visual=subplot('position',[.03 .66 .26 .29]);
    Ax_visual_bar=subplot('position',[.03 .63 .26 .02]);
    
    startx = uicontrol(f,'style','pushbutton','units','normalized','position',[0.4690    0.58    0.0520    0.05],'string','Start','fontsize',15,'callback',@fcn_start);      
    queryModeH=uicontrol(f,'Style','checkbox','String','Query mode','Value',1,'units','normalized','Position',[0.93    0.62    0.05   0.025],'HorizontalAlignment','left','callback',@fcn_queryMode);
    queryMode=1;
    autoNextH=uicontrol(f,'Style','checkbox','String','Auto-next','Value',1,'units','normalized','Position',[0.93    0.64    0.05   0.025],'HorizontalAlignment','left','callback',@fcn_autoNext);
    autoNext=1;
    tm=uicontrol(f,'style','pushbutton','units','normalized','position',[0.93    0.67    0.05   0.025],'string','Template match','callback',@fcn_tm);
    trainRF=uicontrol(f,'style','pushbutton','units','normalized','position',[0.93    0.59    0.05   0.025],'string','Next query','callback',@fcn_trainRF);
    SeqMode='Sequential off';
    seqModePop=uicontrol('Style','popup','String',{'Sequential off','Seizure','LPD','GPD','LRDA','GRDA','Other'},'units','normalized','Position',[0.93    0.56    0.05   0.025],'value',1,'Callback',@fcn_seqModePop);
    ModMode='Model off';
    modModePop=uicontrol('Style','popup','String',{'Model off','Seizure','LPD','GPD','LRDA','GRDA','Other'},'units','normalized','Position',[0.93    0.53    0.05   0.025],'value',1,'Callback',@fcn_modModePop);
    done=uicontrol(f,'style','pushbutton','units','normalized','position',[0.93    0.45    0.05   0.025],'string','Done','callback',@fcn_done);
    bg=uibuttongroup(f,'Title','EEG patterns','units','normalized','Position',[.93 .7 .06 .25],'SelectionChangedFcn',@fcn_decision);
    patterns={'Seizure','LPD','GPD','LRDA','GRDA','Other'};   
    nP=length(patterns);
    r=cell(nP,1);
    ro=uicontrol(bg,'Style','radiobutton','String','None','units','normalized','Position',[.1 1 .8 .1],'HandleVisibility','off');
    for ii=1:nP
        r{ii}=uicontrol(bg,'Style','radiobutton','String',['[',shortcuts{ii},'] ',patterns{ii}],'units','normalized','Position',[.1 (-ii+nP+1)/(1+nP) .8 .08],'HandleVisibility','off');
    end
    firstSeg=uicontrol(f,'style','pushbutton','units','normalized','position',[.030 .58 .015 .025],'string','1st','callback',@fcn_1stSeg);
    prevSeg=uicontrol(f,'style','pushbutton','units','normalized','position',[.045 .58 .03 .025],'string','Previous','callback',@fcn_prevSeg);
    nextSeg=uicontrol(f,'style','pushbutton','units','normalized','position',[.075 .58 .033 .025],'string','Next epoch','callback',@fcn_nextSeg);
    firstSZ=uicontrol(f,'style','pushbutton','units','normalized','position',[.14 .58 .015 .025],'string','1st','callback',@fcn_1stSZ);
    prevSZ=uicontrol(f,'style','pushbutton','units','normalized','position',[.155 .58 .03 .025],'string','Previous','callback',@fcn_prevSZ);
    nextSZ=uicontrol(f,'style','pushbutton','units','normalized','position',[.185 .58 .03 .025],'string','Next SZ','callback',@fcn_nextSZ);
    winSize=10;
    winSizeH=uicontrol('Style','popup','String',{'60 min','30 min','10 min','2 min'},'units','normalized','Position',[.26 .58 .03 .025],'value',3,'Callback',@fcn_winSize);
    colorMode='By pattern';
    colorPop=uicontrol('Style','popup','String',{'By pattern'},'units','normalized','Position',[.23 .95 .06 .025],'value',1,'Callback',@fcn_colorPop);
    montage='L-Bipolar';
    Montage=uicontrol('Style','popup','String',{'L-Bipolar','Average','Monopolar'},'units','normalized','Position',[.3 .95 0.045 0.0250],'value',1,'Callback',@fcn_montageFlag);
    
    set(0,'CurrentFigure',f1)
    Ax_spec1={subplot('position',[.03 .745 .26 .205]);
                  subplot('position',[.03 .530 .26 .205]);
                  subplot('position',[.03 .315 .26 .205]);
                  subplot('position',[.03 .100 .26 .205]);
                  subplot('position',[.03 .950 .26 .010])};
    Ax_EEG1=subplot('position',[.32 .1 .60 .85]);
    bg1=uibuttongroup(f1,'Title','EEG patterns','units','normalized','Position',[.93 .7 .06 .25],'SelectionChangedFcn',@fcn_tmdecision);  
    r1=cell(nP,1);
    ro1=uicontrol(bg1,'Style','radiobutton','String','None','units','normalized','Position',[.1 1 .8 .1],'HandleVisibility','off');
    for ii=1:nP
        r1{ii}=uicontrol(bg1,'Style','radiobutton','String',['[',shortcuts{ii},'] ',patterns{ii}],'units','normalized','Position',[.1 (-ii+nP+1)/(1+nP) .8 .08],'HandleVisibility','off');
    end
    goBack=uicontrol(f1,'style','pushbutton','units','normalized','position',[0.93    0.67    0.05   0.025],'string','Go back','callback',@fcn_goBack);
    
    set(Ax_EEG,'Visible','off');
    for ii=1:length(Ax_spec)
        set(Ax_spec{ii},'Visible','off');
    end
    set(Ax_visual,'Visible','off');
    set(Ax_visual_bar,'Visible','off');
    set(bg,'Visible','off');
    set(trainRF,'Visible','off');
    set(done,'Visible','off');
    set(colorPop,'Visible','off');
    set(winSizeH,'Visible','off');
    set(Montage,'Visible','off');
    set(autoNextH,'Visible','off');
    set(queryModeH,'Visible','off');
    set(firstSeg,'Visible','off');
    set(prevSeg,'Visible','off');
    set(nextSeg,'Visible','off');
    set(firstSZ,'Visible','off');
    set(prevSZ,'Visible','off');
    set(nextSZ,'Visible','off');
    set(tm,'Visible','off');
    set(seqModePop,'Visible','off');
    set(modModePop,'Visible','off');
    
    iter=[];
    N=10;
    jj=1;
    LUTtopN_current=cell(N,6);
    zScale=1/150;
    idxReal=cell(0,8);
    indTopN=[];
    LUT_previous=[];
    doneFlag=[];
    tmFlag=0;
    LUTtopN=[];
    seenFlag=[];
    Xvisual=[];
    bow_vec=[];
    indSeen_=[];
    indUnseen_=[];
    nSamples=[];
    lonerId=0;
    jjj=[];
    jjj1=[];
    y_tm='';
    Y_model=[];
    cluster_lut=[];
    idx_last=[];
    Fs=200;
    col=[-10 25];
    spatialRegs={'LL','RL','LP','RP'};
    channel_withspace_bipolar={'Fp1-F7' 'F7-T3' 'T3-T5' 'T5-O1' '' 'Fp2-F8' 'F8-T4' 'T4-T6' 'T6-O2' '' 'Fp1-F3' 'F3-C3' 'C3-P3' 'P3-O1' '' 'Fp2-F4' 'F4-C4' 'C4-P4' 'P4-O2' '' 'Fz-Cz'  'Cz-Pz' '' 'EKG'};
    channel_withspace_average={'Fp1-avg' 'F3-avg' 'C3-avg' 'P3-avg' 'F7-avg' 'T3-avg' 'T5-avg' 'O1-avg' '' 'Fz-avg' 'Cz-avg' 'Pz-avg' '' 'Fp2-avg' 'F4-avg' 'C4-avg' 'P4-avg' 'F8-avg' 'T4-avg' 'T6-avg' 'O2-avg' '' 'EKG'};
    channel_withspace_monopolar={'Fp1' 'F3' 'C3' 'P3' 'F7' 'T3' 'T5' 'O1' '' 'Fz' 'Cz' 'Pz' '' 'Fp2' 'F4' 'C4' 'P4' 'F8' 'T4' 'T6' 'O2' '' 'EKG'};
    
    uiwait(f);
    while true
        set(done,'Enable','on');
        set(colorPop,'Enable','on');
        set(Montage,'Enable','on');
        set(autoNextH,'Enable','on');
        set(queryModeH,'Enable','on');
        set(firstSeg,'Enable','on');
        set(prevSeg,'Enable','on');
        set(nextSeg,'Enable','on');
        set(firstSZ,'Enable','on');
        set(prevSZ,'Enable','on');
        set(nextSZ,'Enable','on');
        set(trainRF,'Enable','on');
        set(tm,'Enable','on');
        set(seqModePop,'Enable','on')
        set(modModePop,'Enable','on')        
        for ii=1:nP
            set(r{ii},'Enable','on')
        end
        set(ro,'Enable','on');
        set(0,'CurrentFigure',f)
        
        if strcmp(SeqMode,'Sequential off')==0||strcmp(ModMode,'Model off')==0  
            y=LUT_previous{jjj,4};
            y_hat=[];
            iEpoch=LUT_previous{jjj,1};
            ii=LUT_previous{jjj,2};
        else
            if queryMode==1                    
                y=LUTtopN{jj,4};
                y_hat=LUTtopN_current{jj,4};
                iEpoch=LUTtopN{jj,1};
                ii=LUTtopN{jj,2};
                jjj=indTopN(jj);
            else  
                y=LUT_previous{jjj,4};
                y_hat=[];
                iEpoch=LUT_previous{jjj,1};
                ii=LUT_previous{jjj,2};
            end
        end     
        
        tmp=load([dataDir,iEpoch,'_',num2str(ii),'.mat']);
        seg=tmp.SEG;sdata=tmp.Sdata;sfreqs=tmp.sfreqs;timeStampo=datetime('2001-01-01 00:00:00');
        eeg=seg(1:19,:);
        eeg(eeg>300)=300;eeg(eeg<-300)=-300;
        ekg=seg(20,:);
        
        w=14;tc=ii*2-1;tto=tc*Fs-(w/2)*Fs+1;tt1=tc*Fs+(w/2)*Fs;tt=tto:tt1;
        gap=NaN(1,size(eeg,2));
        switch montage
            case 'L-Bipolar'
                seg=fcn_bipolar(eeg);
                seg_disp=[seg(1:4,:);gap;seg(5:8,:);gap;seg(9:12,:);gap;seg(13:16,:);gap;seg(17:18,:);gap;ekg];
                channel_withspace=channel_withspace_bipolar;
            case 'Average' 
                seg=eeg-repmat(mean(eeg,1),size(eeg,1),1);
                seg_disp=[seg(1:8,:);gap;seg(9:11,:);gap;seg(12:19,:);gap;ekg];
                channel_withspace=channel_withspace_average;
            case 'Monopolar'
                seg=eeg;
                seg_disp=[seg(1:8,:);gap;seg(9:11,:);gap;seg(12:19,:);gap;ekg];
                channel_withspace=channel_withspace_monopolar;
        end
        M=size(seg_disp,1);DCoff=repmat(flipud((1:M)'),1,size(seg_disp,2));
        ekg_=seg_disp(end,:);ekg_=(ekg_-mean(ekg_))/(eps+std(ekg_));
        timeStamps=datestr(timeStampo+seconds(round(tto/Fs:2:(tt1+1)/Fs)),'hh:MM:ss');
        
        set(f,'CurrentAxes',Ax_EEG);cla(Ax_EEG)
        hold(Ax_EEG,'on')
            for ii=1:round((tt1-tto+1)/Fs)
                ta=tto+Fs*(ii-1);
                line([ta ta],[0 M+1],'linestyle','--','color',[.5 .5 .5])
            end
            plot(Ax_EEG,tt,zScale*seg_disp(1:end-1,:)+DCoff(1:end-1,:),'k');
            plot(Ax_EEG,tt,-.2*ekg_+DCoff(end,:),'m');           
            set(Ax_EEG,'ytick',1:M,'yticklabel',fliplr(channel_withspace),'box','on','ylim',[0 M+1],'xlim',[tto tt1+1],'xtick',round(tt(1):2*Fs:(tt(end)+Fs)),'xticklabel',timeStamps)
            
            dt=tt1-tto+1;a=round(dt*4/5);
            xa1=tto+[a a+Fs-1];ya1=[3 3];xa2=tto+[a a];ya2=ya1+[0 100*zScale];
            text(Ax_EEG,xa1(1)-.7*a/10,mean(ya2),'100\muV','Color','r','FontSize',12);text(Ax_EEG,mean(xa1)-0.3*a/10,2.5,'1 sec','Color','r','FontSize',12);
            line(Ax_EEG,xa1,ya1,'LineWidth',2,'Color','r');line(Ax_EEG,xa2,ya2,'LineWidth',2,'Color','r');            
            
            if queryMode==1
                if doneFlag(jj)==1
                    set(bg,'SelectedObject',r{ismember(patterns,y_hat)})
                    title(['[',seenFlag{jj},'] ',' Pseudo [\color{red}',y,'\color{black}] vs.',' Expert [\color{blue}',y_hat,'\color{black}]'],'FontSize',20)
                else  
                    set(bg,'SelectedObject',ro)
                    title(['[',seenFlag{jj},'] ',' Pseudo [\color{red}',y,'\color{black}] vs.',' Expert [\color{blue}','\color{black}]'],'FontSize',20)
                end
            else 
                tmp_=jjj;
                if ~isempty(find(indTopN==tmp_,1))          
                    jj_=find(indTopN==tmp_);
                    if doneFlag(jj_)==1
                        set(bg,'SelectedObject',r{ismember(patterns,LUTtopN_current{jj_,4})})
                        title(['[',seenFlag{jj_},'] ',' Pseudo [\color{red}',LUTtopN_current{jj_,3},'\color{black}] vs.',' Expert [\color{blue}',LUTtopN_current{jj_,4},'\color{black}]'],'FontSize',20)
                    else  
                        set(bg,'SelectedObject',ro)
                        title(['[',seenFlag{jj_},'] ',' Pseudo [\color{red}',LUTtopN_current{jj_,3},'\color{black}] vs.',' Expert [\color{blue}','\color{black}]'],'FontSize',20)
                    end
                else
                    set(bg,'SelectedObject',ro)
                     title(['Pseudo [\color{red}',y,'\color{black}] vs.',' Expert [\color{blue}','\color{black}]'],'FontSize',20)
                end 
            end
        hold(Ax_EEG,'off')
        
        winSize=10;stepSize=2;
        to=(tc-1)-(winSize*60/2)+1;t1=(tc-1)+(winSize*60/2);
        S_x=to:stepSize:t1;S_y=sfreqs;
        for ii=1:size(sdata,1)
            set(f,'CurrentAxes',Ax_spec{ii});cla(Ax_spec{ii});
            colormap jet;axis(Ax_spec{ii},'xy'); 
            spec=sdata{ii,1};
            hold(Ax_spec{ii},'on')
                imagesc(Ax_spec{ii},S_x,S_y,pow2db(spec),col);
                plot([tc  tc],[S_y(1) S_y(end)],'k--','linewidth',1);
            hold(Ax_spec{ii},'off')
            ylabel(Ax_spec{ii},'Freq (Hz)')
            ss=get(Ax_spec{ii},'yticklabel');ss{end}=spatialRegs{ii};
            set(Ax_spec{ii},'xtick',[],'yticklabel',ss,'ylim',[S_y(1) S_y(end)],'xlim',[tc-winSize/2*60,tc+winSize/2*60],'box','on') 
        end
        
        tt=timeStampo+seconds(to:2*60:(t1+1));
        hold(Ax_spec{4},'on')  
            set(Ax_spec{4},'xtick',to:2*60:(t1+1),'xticklabel',datestr(tt,'hh:MM:ss'),'xlim',[tc-winSize/2*60,tc+winSize/2*60]);
        hold(Ax_spec{4},'off') 
        
        set(f,'CurrentAxes',Ax_spec{end});cla(Ax_spec{end})
        plot(Ax_spec{end},tc,0,'rv','markersize',8,'MarkerFaceColor','r');
        axis(Ax_spec{end},'off');
        set(Ax_spec{end},'xtick',[],'xlim',get(Ax_spec{end-1},'xlim'),'ylim',[0 .5]) 
        
        set(f,'CurrentAxes',Ax_visual_bar); 
        colors=flipud(jet(7));colors=colors(1:6,:);colors=reshape(colors,1,size(colors ,1),3);
        imagesc(Ax_visual_bar,1:6,1,colors)
        set(Ax_visual_bar,'xtick',1:6,'xticklabel',patterns,'ytick',[])
        
        set(f,'CurrentAxes',Ax_visual);cla(Ax_visual);
        Vxy=Xvisual;
        lut=LUT_previous;lut(indTopN(doneFlag==1),4)=LUTtopN_current(doneFlag==1,4);
        colors=fcn_pattern2color(lut,colorMode);
        hold(Ax_visual,'on');
            ss=scatter(Ax_visual,Vxy(:,1),Vxy(:,2),20,colors(:,:),'filled');alpha(ss,.2)
            scatter(Ax_visual,Vxy(indSeen_,1),Vxy(indSeen_,2),20,colors(indSeen_,:),'filled');
            if strcmp(SeqMode,'Sequential off')==0 ||strcmp(ModMode,'Model off')==0
                plot(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),'mx','markersize',20,'linewidth',2);plot(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),'mo','markersize',20,'linewidth',2);
                scatter(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),20,colors(jjj,:),'filled');
                
                id_done=indTopN(doneFlag==1);
                scatter(Ax_visual,Vxy(id_done,1),Vxy(id_done,2),20,colors(id_done,:),'v','filled','MarkerEdgeColor','k');
                text(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),['    pseudo: ',y,newline,'    expert : ',y_hat]);
            else
                if queryMode==1 
                    plot(Ax_visual,Vxy(indTopN(jj),1),Vxy(indTopN(jj),2),'mx','markersize',20,'linewidth',2);plot(Ax_visual,Vxy(indTopN(jj),1),Vxy(indTopN(jj),2),'mo','markersize',20,'linewidth',2);
                    scatter(Ax_visual,Vxy(indTopN(jj),1),Vxy(indTopN(jj),2),20,colors(indTopN(jj),:),'filled');
                   
                    id_done=indTopN(doneFlag==1);
                    scatter(Ax_visual,Vxy(id_done,1),Vxy(id_done,2),20,colors(id_done,:),'v','filled','MarkerEdgeColor','k');
                    text(Ax_visual,Vxy(indTopN(jj),1),Vxy(indTopN(jj),2),['    pseudo: ',y,newline,'    expert : ',y_hat]);
                else
                    plot(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),'mx','markersize',20,'linewidth',2);plot(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),'mo','markersize',20,'linewidth',2);
                    scatter(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),20,colors(jjj,:),'filled');
                    
                    id_done=indTopN(doneFlag==1);% ind in LUT %
                    scatter(Ax_visual,Vxy(id_done,1),Vxy(id_done,2),20,colors(id_done,:),'v','filled','MarkerEdgeColor','k');
                    text(Ax_visual,Vxy(jjj,1),Vxy(jjj,2),['    pseudo: ',y,newline,'    expert : ',y_hat]);
                end
            end
            axis(Ax_visual,'off');axis(Ax_visual,'equal');set(get(Ax_visual,'title'),'string',['iter',num2str(iter)])
        hold(Ax_visual,'off')
        uiwait(f);
    end

    function fcn_start(varargin) 
        set(startx,'Visible','off','Enable','off');
        loading=uicontrol(f,'style','text','units','normalized','position',[0.4690 0.5800 0.0520 0.0500],'string','Loading...','FontSize',15);
        drawnow;    
        if exist(targetDir,'dir')==7
            al_files=dir([targetDir,'ActiveLearning_LUT_round*']);
            iter=size(al_files,1);
            if iter==0
                tmp=load('.\Data\PaCMAP\pacmap_output.mat');
                Y_model=tmp.Y_model;
                Xvisual=tmp.Vxy;               
                tmp=load('.\Data\BoW.mat');
                bow_vec=tmp.bow_vec;
                tmp=load('.\Data\GUI_LUT.mat');
                LUT=tmp.LUT_all;
                LUT=LUT(:,[1 2 4 4]);
                save([targetDir,'ActiveLearning_LUT_round0.mat'],'LUT','Y_model','Xvisual','bow_vec')
                iter=1;
            end  
        else  
            mkdir(targetDir)
            tmp=load('.\Data\PaCMAP\pacmap_output.mat');
            Y_model=tmp.Y_model;
            Xvisual=tmp.Vxy;
            tmp=load('.\Data\BoW.mat');
            bow_vec=tmp.bow_vec;            
            tmp=load('.\Data\GUI_LUT.mat');
            LUT=tmp.LUT_all;
            LUT=LUT(:,[1 2 4 4]); 
            save([targetDir,'ActiveLearning_LUT_round0.mat'],'LUT','Y_model','Xvisual','bow_vec')
            iter=1;
        end
        if exist([targetDir,'ActiveLearning_LUT_final.mat'],'file')==2
            tmp_previous=load([targetDir,'ActiveLearning_LUT_final.mat']);
            LUT_previous=tmp_previous.LUT;             
            LUT_previous=[LUT_previous,LUT_previous(:,end)];
            indSeen_=tmp_previous.idxSeen;
            Y_model=tmp_previous.Y_model;
            Xvisual=tmp_previous.Xvisual;
            bow_vec=tmp_previous.bow_vec;
            idxReal=tmp_previous.idxReal;
        else
            tmp_previous=load([targetDir,'ActiveLearning_LUT_round',num2str(iter-1),'.mat']);
            LUT_previous=tmp_previous.LUT;
            Y_model=tmp_previous.Y_model;
            Xvisual=tmp_previous.Xvisual;
            bow_vec=tmp_previous.bow_vec;
            if iter==1
                indSeen_=[];
                idxReal=cell(0,8);
            else 
                indSeen_=tmp_previous.indSeen_;
                idxReal=tmp_previous.idxReal;
            end
        end 
        nSamples=size(LUT_previous,1);
        indUnseen_=setdiff((1:nSamples)',indSeen_);
        if isempty(indUnseen_)
            choice=questdlg('Wow! Seems you have labeled all the segments!','Warning','Ok,I am done.','Ok,I am done.');
            switch choice
                case 'Ok,I am done.'
                    fcn_done();
            end
            close all
            delete(f)           
        else
            fcn_dataQuery_model();
            N=length(indTopN);
            LUTtopN_current=cell(N,4);
            LUTtopN_current(:,1:3)=LUTtopN(:,[1:2,4]);
            doneFlag=NaN(N,1);
            jj=1;
            jjj=indTopN(jj);
            delete(loading)
            set(Ax_EEG,'Visible','on');
            for iAx=1:length(Ax_spec)
                set(Ax_spec{iAx},'Visible','on');
            end
            set(Ax_visual,'Visible','on');
            set(Ax_visual_bar,'Visible','on');
            set(bg,'Visible','on');
            set(trainRF,'Visible','on','Enable','off');
            set(done,'Visible','on');
            set(colorPop,'Visible','on');
            set(Montage,'Visible','on');
            set(autoNextH,'Visible','on');
            set(queryModeH,'Visible','on');
            set(firstSeg,'Visible','on');
            set(prevSeg,'Visible','on');
            set(nextSeg,'Visible','on');
            set(firstSZ,'Visible','on');
            set(prevSZ,'Visible','on');
            set(nextSZ,'Visible','on');
            set(tm,'Visible','on')
            set(seqModePop,'Visible','on')
            set(modModePop,'Visible','on') 
            uiresume(f);
        end
    end

    function dataBipolar=fcn_bipolar(data)
        dataBipolar( 1,:)=data( 1,:)-data( 5,:);% Fp1-F7
        dataBipolar( 2,:)=data( 5,:)-data( 6,:);% F7-T3
        dataBipolar( 3,:)=data( 6,:)-data( 7,:);% T3-T5
        dataBipolar( 4,:)=data( 7,:)-data( 8,:);% T5-O1

        dataBipolar( 5,:)=data(12,:)-data(16,:);% Fp2-F8
        dataBipolar( 6,:)=data(16,:)-data(17,:);% F8-T4
        dataBipolar( 7,:)=data(17,:)-data(18,:);% T4-T6
        dataBipolar( 8,:)=data(18,:)-data(19,:);% T6-O2
        
        dataBipolar( 9,:)=data( 1,:)-data( 2,:);% Fp1-F3
        dataBipolar(10,:)=data( 2,:)-data( 3,:);% F3-C3
        dataBipolar(11,:)=data( 3,:)-data( 4,:);% C3-P3
        dataBipolar(12,:)=data( 4,:)-data( 8,:);% P3-O1

        dataBipolar(13,:)=data(12,:)-data(13,:);% Fp2-F4
        dataBipolar(14,:)=data(13,:)-data(14,:);% F4-C4
        dataBipolar(15,:)=data(14,:)-data(15,:);% C4-P4
        dataBipolar(16,:)=data(15,:)-data(19,:);% P4-O2

        dataBipolar(17,:)=data( 9,:)-data(10,:);% Fz-Cz
        dataBipolar(18,:)=data(10,:)-data(11,:);% Cz-Pz
    end

    function keys_Callback(f,varargin)
        k=get(f,'CurrentKey');
        switch k
            case 'rightarrow'             
                if tmFlag==1 
                    tmFlag=0;set(f1,'Visible','off');
                    autoNext=1;set(autoNextH,'value',autoNext);    
                else
                    if queryMode==1 
                        jj=jj+1;
                        if jj>length(indTopN)
                            jj=1;
                        end 
                    else
                        if strcmp(SeqMode,'Sequential off')==0 
                            fcn_sequentialReviewByPattern(SeqMode);
                            
                        elseif strcmp(ModMode,'Model off')==0
                            fcn_modelRecommendByPattern(ModMode);
                        else
                            jjj=jjj+1;
                            if jjj>size(LUT_previous,1)
                                jjj=1;
                            end
                        end
                    end
                end               
            case 'leftarrow'
                if tmFlag==1
                    tmFlag=0;
                    autoNext=1;set(autoNextH,'value',autoNext);
                    set(f1,'Visible','off');
                else
                    if queryMode==1 
                        jj=jj-1;
                        if jj<1
                            jj=length(indTopN);
                        end
                    else
                        if strcmp(SeqMode,'Sequential off')==0 
                            if ~isempty(idx_last)
                                jjj=idx_last;
                                idx_last=[];
                            else
                                fcn_sequentialReviewByPattern(SeqMode);
                            end                        
                        elseif strcmp(ModMode,'Model off')==0
                            if ~isempty(idx_last)
                                jjj=idx_last;
                                idx_last=[];
                            else
                                fcn_modelRecommendByPattern(ModMode);
                            end                          
                        else
                            jjj=jjj-1;
                            if jjj<1
                                jjj=nSamples;
                            end
                        end
                    end
                end
            case 'uparrow'
                zScale=zScale*1.5;          
            case 'downarrow'
                zScale=zScale/1.5;   
            case {'1','numpad1'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{1})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{1})
                    fcn_decisionkeys(1);
                end
            case {'2','numpad2'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{2})
                    fcn_tmdecisionkeys(2);
                else
                    set(bg,'SelectedObject',r{2})
                    fcn_decisionkeys(2);
                end             
            case {'3','numpad3'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{3})
                    fcn_tmdecisionkeys(3);
                else
                    set(bg,'SelectedObject',r{3})
                    fcn_decisionkeys(3);
                end         
            case {'4','numpad4'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{4})
                    fcn_tmdecisionkeys(4);
                else
                    set(bg,'SelectedObject',r{4})
                    fcn_decisionkeys(4);
                end        
            case {'5','numpad5'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{5})
                    fcn_tmdecisionkeys(5);
                else
                    set(bg,'SelectedObject',r{5})
                    fcn_decisionkeys(5);
                end     
            case {'6','numpad6'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{6})
                    fcn_tmdecisionkeys(6);
                else
                    set(bg,'SelectedObject',r{6})
                    fcn_decisionkeys(6);
                end  
            case {'space'}
                if tmFlag==1
                    tmFlag=0;
                    autoNext=1;set(autoNextH,'value',autoNext);
                    set(f1,'Visible','off');
                else
                    fcn_nextSeg();
                end
            case {'escape'}
                if tmFlag==1
                    fcn_goBack(); 
                end               
            otherwise
                disp(['Dude! Invalid Key ',k])
        end
        uiresume(f);
    end

    function fcn_decisionkeys(ind_)   
        idx_last=jjj;
        p=patterns{ind_};
        idxReal=[idxReal;{jjj,p,datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF'),autoNext,queryMode,tmFlag,SeqMode,ModMode}];
        tmp_=jjj;
        if isempty(find(indTopN==tmp_,1))        
            lonerId=lonerId+1;
            indTopN=[indTopN;tmp_];
            LUTtopN_current=[LUTtopN_current;[LUT_previous(tmp_,[1:2 4]),{[]}]];
            LUTtopN=[LUTtopN;LUT_previous(tmp_,:)];
            doneFlag=[doneFlag;NaN];
            seenFlag=[seenFlag;['Free point ',num2str(lonerId)]];
            jj_=length(indTopN);
        else
            jj_=find(indTopN==tmp_);% ind in indTopN_ %
        end
        LUTtopN_current{jj_,4}=p;
        doneFlag(jj_)=1;
        if strcmp(SeqMode,'Sequential off')==0
            fcn_sequentialReviewByPattern(SeqMode);        
        else
            if strcmp(ModMode,'Model off')==0
                fcn_modelRecommendByPattern(ModMode);         
            else 
                if tmFlag~=0
                    fcn_tm();
                    %uiresume(f);
                else
                    if autoNext&&queryMode 
                        fcn_trainRF();
                    else
                        fcn_nextCPC();
                    end
                end 
            end
        end
    end

    function fcn_tmdecisionkeys(ind_)
        p=patterns{ind_};
        idxReal=[idxReal;{jjj1,p,datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF'),autoNext,queryMode,tmFlag,SeqMode,ModMode}];    
        tmp_=jjj1;     
        if isempty(find(indTopN==tmp_,1))       
            lonerId=lonerId+1;
            indTopN=[indTopN;tmp_];
            LUTtopN_current=[LUTtopN_current;[LUT_previous(tmp_,[1:2 4]),{[]}]];
            LUTtopN=[LUTtopN;LUT_previous(tmp_,:)];
            doneFlag=[doneFlag;NaN];
            seenFlag=[seenFlag;['Free point ',num2str(lonerId)]];
            jj_=length(indTopN);
        else
            jj_=find(indTopN==tmp_); 
        end
        LUTtopN_current{jj_,4}=p;
        doneFlag(jj_)=1;
        if strcmp(y_tm,p)
            fcn_tm();
        else
            tmFlag=0;
            autoNext=1;set(autoNextH,'value',autoNext);
            set(f1,'Visible','off');
            uiresume(f);
        end       
    end

    function fcn_decision(source,event)
        for ir=1:nP
            set(r{ir},'Enable','off');
        end
        set(ro,'Enable','off');
        drawnow;    
        idx_last=jjj;   
        p=event.NewValue.String;ip=regexp(p,']');p=p(ip+2:end);
        idxReal=[idxReal;{jjj,p,datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF'),autoNext,queryMode,tmFlag,SeqMode,ModMode}];      
        tmp_=jjj;
        if isempty(find(indTopN==tmp_,1))         
            lonerId=lonerId+1;
            indTopN=[indTopN;tmp_];
            LUTtopN_current=[LUTtopN_current;[LUT_previous(tmp_,[1:2 4]),{[]}]];
            LUTtopN=[LUTtopN;LUT_previous(tmp_,:)];
            doneFlag=[doneFlag;NaN];
            seenFlag=[seenFlag;['Free point ',num2str(lonerId)]];
            jj_=length(indTopN);
        else
            jj_=find(indTopN==tmp_);
        end
        LUTtopN_current{jj_,4}=p;
        doneFlag(jj_)=1;
        if strcmp(SeqMode,'Sequential off')==0 
            fcn_sequentialReviewByPattern(SeqMode);        
        else
            if strcmp(ModMode,'Model off')==0
                fcn_modelRecommendByPattern(ModMode);          
            else
                if tmFlag~=0 
                    fcn_tm();
                    %uiresume(f);
                else
                    if autoNext&&queryMode
                        fcn_trainRF();
                    else
                        fcn_nextCPC();
                    end
                end
            end
        end
    end

    function fcn_tmdecision(source,event)
        for ir=1:nP
            set(r1{ir},'Enable','off');
        end
        set(ro1,'Enable','off');
        drawnow;        
        tmFlag=1;set(f1,'Visible','on');    
        autoNext=0;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);        
        p=event.NewValue.String;
        ip=regexp(p,']');
        p=p(ip+2:end);
        idxReal=[idxReal;{jjj1,p,datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF'),autoNext,queryMode,tmFlag,SeqMode,ModMode}];
        tmp_=jjj1;      
        if isempty(find(indTopN==tmp_,1))        
            lonerId=lonerId+1;
            indTopN=[indTopN;tmp_];
            LUTtopN_current=[LUTtopN_current;[LUT_previous(tmp_,[1:2 4]),{[]}]];
            LUTtopN=[LUTtopN;LUT_previous(tmp_,:)];
            doneFlag=[doneFlag;NaN];
            seenFlag=[seenFlag;['Free point ',num2str(lonerId)]];
            jj_=length(indTopN);
        else
            jj_=find(indTopN==tmp_);
        end
        LUTtopN_current{jj_,4}=p;
        doneFlag(jj_)=1;
        if strcmp(y_tm,p)
            fcn_tm();
        else
            autoNext=1;set(autoNextH,'value',autoNext);          
            queryMode=0;set(queryModeH,'value',queryMode);
            tmFlag=0;set(f1,'Visible','off');
            SeqMode='Sequential off';set(seqModePop,'value',1);
            ModMode='Model off';set(modModePop,'value',1);  
            uiresume(f);
        end 
    end

    function fcn_goBack(varargin)
        set(goBack,'Enable','off');
        tmFlag=0;
        autoNext=1;set(autoNextH,'value',autoNext);
        set(f1,'Visible','off');
        uiresume(f);
    end

    function fcn_nextCPC(~,~)
        jjj=jjj+1;
        if jjj>size(LUT_previous,1)
            jjj=1;
        end
        uiresume(f);
    end

    function fcn_trainRF(varargin)        
        set(trainRF,'Enable','off');
        drawnow
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=1;set(queryModeH,'value',queryMode);       
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
       
        if ~isempty(find(isnan(doneFlag)==1,1)) 
            jj=find(isnan(doneFlag)==1,1);
            uiresume(f);
        
        else
            choice=questdlg('Update the MAP?','Warning','Yes','Yes');
            switch choice
                case 'Yes' 
                    fcn_spreadingCallback();
                    iter=iter+1;
                    jj=1;
                    lonerId=0; 
                    set(Ax_EEG,'Visible','off');cla(Ax_EEG)
                    for iAx=1:length(Ax_spec)
                        set(Ax_spec{iAx},'Visible','off');
                        cla(Ax_spec{iAx})
                    end
                    set(Ax_visual,'Visible','off');cla(Ax_visual)
                    set(Ax_visual_bar,'Visible','off');cla(Ax_visual_bar)
                    set(get(Ax_visual,'title'),'string','')
                    set(bg,'Visible','off');
                    set(trainRF,'Visible','off');
                    set(done,'Visible','off');
                    set(colorPop,'Visible','off');
                    set(Montage,'Visible','off');
                    set(autoNextH,'Visible','off');
                    set(queryModeH,'Visible','off');
                    set(firstSeg,'Visible','off');
                    set(prevSeg,'Visible','off');
                    set(nextSeg,'Visible','off');
                    set(firstSZ,'Visible','off');
                    set(prevSZ,'Visible','off');
                    set(nextSZ,'Visible','off');
                    set(tm,'Visible','off');
                    set(seqModePop,'Visible','off');
                    set(modModePop,'Visible','off');
                    
                    drawnow
                    fcn_start();
            end
        end       
    end

    function fcn_dataQuery_model(~,~)
        y_expert=fcn_label2number([],categorical(LUT_previous(:,4))); 
        [yp_model,y_model]=max(Y_model(:,[2:6,1]),[],2);
        con_unseen=yp_model(indUnseen_);
        scr_unseen=y_model(indUnseen_);
        tru_unseen=y_expert(indUnseen_);
        Xunseen_=Xvisual(indUnseen_ ,:);        
        scr_seen=y_model(indSeen_);
        tru_seen=y_expert(indSeen_);
        Xseen_=Xvisual(indSeen_ ,:);        
        K=4;
        indTopN=[];
        seenFlag=cell(0,1);
        if iter==1||isempty(indSeen_) 
            cluster_lut=cell(6,3); 
            for p=1:6  
                idx_p=find(scr_unseen==p);
                if length(idx_p)<=K  
                    con_p=con_unseen(idx_p);
                    [~,idx1]=sort(con_p,'descend'); 
                    indTopN_p=indUnseen_(idx_p(idx1)); 
                    indTopN=[indTopN;indTopN_p]; 
                    cluster_lut(p,:)={p,indTopN_p,[]};                    
                else
                    Xunseen_p=Xunseen_(idx_p,:);
                    [IDX_p,xMedoids_p]=kmedoids(Xunseen_p,K); 
                    con_p=NaN(K,1); 
                    idx_medoids_p=NaN(K,1); 
                    idx_members_p=IDX_p; 
                    for k=1:K
                        x_medoid=xMedoids_p(k,:);
                        idx_k=find(ismember(Xunseen_p,x_medoid,'rows')==1,1);
                        con_p(k)=con_unseen(idx_p(idx_k));
                        idx_medoids_p(k)=indUnseen_(idx_p(idx_k));
                        idx_members_p(IDX_p==k)=idx_medoids_p(k);
                    end
                    cluster_lut(p,:)={p,idx_medoids_p,[indUnseen_(idx_p),idx_members_p]};
                    [~,idx1]=sort(con_p,'descend');
                    indTopN_p=idx_medoids_p(idx1);
                    indTopN=[indTopN;indTopN_p];
                end
                seenFlag_p=cell(length(indTopN_p),1);
                for i=1:length(indTopN_p)
                    seenFlag_p{i}=[patterns{p},'-C',num2str(i)];
                end
                seenFlag=[seenFlag;seenFlag_p];
            end
            for i=1:length(indTopN)
                seenFlag{i}=['Query ',num2str(i),'-',seenFlag{i}];
            end  
        else 
            cluster_lut=[];
            for p=1:6
                idx_unseen_p=find(tru_unseen==p&scr_unseen==p);
                idx_seen_p=find(tru_seen==p  &scr_seen==p);                
                Xunseen_p=Xunseen_(idx_unseen_p,:);
                Xseen_p=Xseen_(idx_seen_p,:);
                nn=length(idx_unseen_p);mm=length(idx_seen_p);
                if nn<=K
                    indTopN_p=indUnseen_(idx_unseen_p);       
                else
                    dd=NaN(nn,mm);
                    for m=1:mm
                        xx=Xseen_p(m,:);
                        dd(:,m)=sum((Xunseen_p-repmat(xx,nn,1)).^2,2);
                    end
                    d=min(dd,[],2);
                    [~,idx]=max(d);
                    indTopN_p=indUnseen_(idx_unseen_p(idx));
                    dd(idx,:)=NaN;Xunseen_p(idx,:)=NaN;
                    while length(indTopN_p)<K&&length(indTopN_p)<nn
                        i=indTopN_p(end);
                        xx=Xvisual(i,:);
                        tmp_dd=sum((Xunseen_p-repmat(xx,nn,1)).^2,2);
                        dd=[dd,tmp_dd];
                        d=min(dd,[],2);
                        [~,idx]=max(d);
                        indTopN_p=[indTopN_p;indUnseen_(idx_unseen_p(idx))];
                        dd(idx,:)=NaN;Xunseen_p(idx,:)=NaN;
                        if length(indTopN_p)>=K||length(indTopN_p)>=nn
                            break
                        end
                    end
                end
                indTopN_p1=indTopN_p;                
                idx_unseen_p=find(tru_unseen~=p&scr_unseen==p);
                if length(idx_unseen_p)<=K
                    indTopN_p=indUnseen_(idx_unseen_p);
                else
                    Xunseen_p=Xunseen_(idx_unseen_p,:);
                    [~,xMedoids_p]=kmedoids(Xunseen_p,K);
                    indTopN_p=NaN(K,1);
                    for k=1:K
                        xx=xMedoids_p(k,:);
                        idx_k=find(ismember(Xunseen_p,xx,'rows')==1,1);
                        indTopN_p(k)=indUnseen_(idx_unseen_p(idx_k));
                    end             
                end
                indTopN_p2=indTopN_p;
                indTopN_p=[indTopN_p1;indTopN_p2];               
                seenFlag_p1=cell(length(indTopN_p1),1);
                for i=1:length(indTopN_p1)
                    seenFlag_p1{i}=[patterns{p},' L',num2str(i)];
                end
                seenFlag_p2=cell(length(indTopN_p2),1);
                for i=1:length(indTopN_p2)
                    seenFlag_p2{i}=[patterns{p},' C',num2str(i)];
                end
                seenFlag_p=[seenFlag_p1;seenFlag_p2];
                con_p=yp_model(indTopN_p);
                [~,idx1]=sort(con_p,'descend');
                indTopN=[indTopN;indTopN_p(idx1)];
                seenFlag=[seenFlag;seenFlag_p(idx1)];
            end            
            for i=1:length(indTopN)
                seenFlag{i}=['Query ',num2str(i),'-',seenFlag{i}];
            end            
        end       
        LUTtopN=LUT_previous(indTopN,:);    
    end

    function fcn_spreadingCallback(~,~)        
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=1;set(queryModeH,'value',queryMode);       
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);       
        indSeen_=unique([indSeen_;indTopN]);
        indUnseen_=setdiff((1:nSamples)',indSeen_);       
        Y_=LUT_previous(:,4);
        Y_(indTopN)=LUTtopN_current(:,4);
        if ~isempty(cluster_lut)&&iter==1
            for p=1:6
                lut_p=cluster_lut{p,3};                
                if ~isempty(lut_p)
                    idx_medoids=cluster_lut{p,2};                   
                    idx_members=lut_p(:,1);
                    y_members=lut_p(:,2);
                    for k=1:length(idx_medoids)
                        idx_k=idx_medoids(k);
                        y_k=Y_{idx_k};
                        idx=find(y_members==idx_k);
                        Y_(idx_members(idx))=repmat({y_k},length(idx),1);
                    end
                end
            end            
        else 
            Xseen_=Xvisual(indSeen_,:);
            Yseen_=Y_(indSeen_);
            Xunseen_=Xvisual(indUnseen_,:);
            N_=length(indUnseen_);M_=length(indSeen_);D_=NaN(N_,M_);
            for m_=1: M_
                x_=Xseen_(m_,:);
                D_(:,m_)=sum((Xunseen_-repmat(x_,N_,1)).^2,2);
            end  
            [~,id_]=min(D_,[],2);
            Yunseen_=Yseen_(id_);
            Y_(indUnseen_)=Yunseen_;
        end        
        Y_new=Y_;
        Y_old=LUT_previous(:,4);
        LUT=LUT_previous;
        LUT(:,3)=Y_old;
        LUT(:,4)=Y_new;
        timeStamp=datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF');
        newLUTname=[targetDir,'ActiveLearning_LUT_round',num2str(iter),'.mat'];
        idxReal_idx=cell2mat(idxReal(:,1));
        [~,idx]=unique(idxReal_idx,'last');
        idxReal=idxReal(idx,:);        
        save(newLUTname,'LUT','Y_model','Xvisual','bow_vec','timeStamp','rater','iter','indSeen_','idxReal');      
        if exist([targetDir,'ActiveLearning_LUT_final.mat'],'file')==2
            idxSeen=indSeen_;
            LUT=LUT(:,[1:2,4]);
            save([targetDir,'ActiveLearning_LUT_final.mat'],'LUT','Xvisual','bow_vec','Y_model','idxSeen','timeStamp','rater','idxReal');
        end
    end

    function c=fcn_pattern2color(LUT,colorMode)
        c=zeros(size(LUT,1),3);
        switch colorMode
            case 'By pattern'
                Cs=flipud(jet(7));            
                identifierList=categorical(LUT(:,end));   
                identifiers=patterns;
                for k=1:length(identifiers)
                    ind=find(identifierList==identifiers{k});
                    c(ind,:)=repmat(Cs(k,:),length(ind),1);
                end
        end
    end

    function clicks_Callback(varargin)
        set(bg,'SelectedObject',ro)     
        xy=get(gca,'CurrentPoint');
        kx=xy(1,1);
        ky=xy(1,2);       
        XLIM=get(Ax_visual,'xlim');
        YLIM=get(Ax_visual,'ylim');
        if  gca==Ax_visual &&(kx>=XLIM(1)&&kx<=XLIM(2)&&ky>=YLIM(1)&&ky<=YLIM(2)) 
            Vxy=Xvisual;
            D=(Vxy(:,1)-kx).^2+(Vxy(:,2)-ky).^2;
            [~,jjj]=min(D);
        end
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
        uiresume(f);       
    end

    function [Ynum,Ystr]=fcn_label2number(Ynum,Ystr)
        if isempty(Ynum)
            Ynum=NaN(length(Ystr),1);
            Ynum(Ystr=='Seizure')=1;
            Ynum(Ystr=='LPD')=2;
            Ynum(Ystr=='GPD')=3;
            Ynum(Ystr=='LRDA')=4;
            Ynum(Ystr=='GRDA')=5;
            Ynum(Ystr=='Other')=6;
        else
            Ystr=cell(length(Ynum),1);
            Ystr(Ynum==1)=repmat({'Seizure'},length(find(Ynum==1)),1);
            Ystr(Ynum==2)=repmat({'LPD'},length(find(Ynum==2)),1);
            Ystr(Ynum==3)=repmat({'GPD'},length(find(Ynum==3)),1);
            Ystr(Ynum==4)=repmat({'LRDA'},length(find(Ynum==4)),1);
            Ystr(Ynum==5)=repmat({'GRDA'},length(find(Ynum==5)),1);
            Ystr(Ynum==6)=repmat({'Other'},length(find(Ynum==6)),1);
        end
    end

    function fcn_colorPop(source,event)
        set(colorPop,'Enable','off');
        drawnow;
        val=source.Value;
        maps=source.String;
        colorMode=maps{val};    
        uiresume(f);
    end

    function fcn_seqModePop(source,event)
        set(seqModePop,'Enable','off');
        drawnow;
        val=source.Value;
        maps=source.String;
        SeqMode=maps{val}; 
        fcn_sequentialReviewByPattern(SeqMode);   
    end

    function fcn_modModePop(source,event)
        set(modModePop,'Enable','off');
        drawnow;
        val=source.Value;
        maps=source.String;
        ModMode=maps{val};  
        fcn_modelRecommendByPattern(ModMode);       
    end

    function fcn_modelRecommendByPattern(p)
        idx_p=find(ismember(patterns,p)==1);
        if isempty(idx_p)
            disp('mod mode off')
            uiresume(f);       
        else
            tmFlag=0;set(f1,'Visible','off');
            autoNext=0;set(autoNextH,'value',autoNext);
            queryMode=0;set(queryModeH,'value',queryMode);
            SeqMode='Sequential off';set(seqModePop,'value',1);
            indSeen_=unique([indSeen_;indTopN(doneFlag==1)]);
            indUnseen_=setdiff((1:nSamples)',indSeen_);
            [yp_model,y_model]=max(Y_model(indUnseen_,[2:6,1]),[],2);
            idx_1=find(y_model==idx_p);
            if isempty(idx_1)
                choice=questdlg('All CP center with this pattern got labeled!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        set(0,'CurrentFigure',f)
                        set(f1,'Visible','off');                     
                        set(modModePop,'value',1);
                        ModMode='Model off';
                        uiresume(f);
                end
            else
                tmp_yp=yp_model(idx_1);
                [~,idx_2]=max(tmp_yp);    
                jjj=indUnseen_(idx_1(idx_2));
                uiresume(f);                
            end
        end
    end

    function fcn_sequentialReviewByPattern(p)
        idx_p=find(ismember(patterns,p)==1);
        if isempty(idx_p)
            disp('seq mode off')
            uiresume(f);           
        else
            tmFlag=0;set(f1,'Visible','off');
            autoNext=0;set(autoNextH,'value',autoNext);
            queryMode=0;set(queryModeH,'value',queryMode);
            ModMode='Model off';set(modModePop,'value',1);
            indSeen_=unique([indSeen_;indTopN(doneFlag==1)]);
            indUnseen_=setdiff((1:nSamples)',indSeen_);  
            y_unseen=fcn_label2number([],categorical(LUT_previous(indUnseen_,4)));
            idx_1=find(y_unseen==idx_p,1);     
            if isempty(idx_1)
                choice=questdlg('All CP center with this pattern got labeled!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        set(0,'CurrentFigure',f)
                        set(f1,'Visible','off');                        
                        set(seqModePop,'value',1);
                        SeqMode='Sequential off';
                        uiresume(f);
                end
            else
                disp(['seq mode on-',p,' jjj=',num2str(jjj)])
                jjj=indUnseen_(idx_1);
                uiresume(f);
            end     
        end
    end

    function fcn_winSize(source,event)
        set(winSizeH,'Enable','off');
        drawnow;
        val=source.Value;
        maps=source.String;       
        winSize=str2double(maps{val}(1:end-4));      
        uiresume(f);
    end

    function fcn_montageFlag(source,event)        
        set(Montage,'Enable','off');
        drawnow;     
        val=source.Value;
        maps=source.String;
        montage=maps{val};        
        uiresume(f);
    end

    function fcn_autoNext(source,event)
        set(autoNextH,'Enable','off');
        drawnow;  
        autoNext=get(autoNextH,'value');
        uiresume(f);
    end
    
    function fcn_queryMode(source,event)
        set(queryModeH,'Enable','off');
        drawnow;        
        queryMode=get(queryModeH,'value');
        if queryMode==1
            autoNext=1;set(autoNextH,'value',autoNext);
            tmFlag=0;set(f1,'Visible','off');
            SeqMode='Sequential off';set(seqModePop,'value',1);
            ModMode='Model off';set(modModePop,'value',1);
        end       
        uiresume(f);
    end

    function fcn_prevSeg(~,~)
        set(prevSeg,'Enable','off');
        drawnow;
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
        jjj_=jjj-(winSize/2*60/2);       
        if jjj_<1
            choice=questdlg('Reach the beginning!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    jjj=1;
            end
        else
            jjj=jjj_;
        end
        uiresume(f);
    end

    function fcn_1stSeg(~,~)
        set(firstSeg,'Enable','off');
        drawnow;
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);    
        jjj=1;
        uiresume(f);
    end

    function fcn_nextSeg(~,~)
        set(nextSeg,'Enable','off');
        drawnow;
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);      
        jjj_=jjj+(winSize/2*60/2);
        if jjj_>size(LUT_previous,1)
            choice=questdlg('Reach the end!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    jjj=size(LUT_previous,1);
            end
        else
            jjj=jjj_;
        end
        uiresume(f);
    end

    function fcn_1stSZ(~,~)
        set(firstSZ,'Enable','off');
        drawnow;
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
        Y_=LUT_previous(:,4);%
        for i_=1:length(indTopN)
            y_=LUTtopN_current{i_,4};
            if ~isempty(y_)
                Y_{indTopN(i_)}=y_;
            end
        end
        idx=find(ismember(Y_,'Seizure')==1);       
        if isempty(idx)
            choice=questdlg('No seizure found!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
        else
            jjj=idx(1);
        end 
 
        uiresume(f);
    end

    function fcn_nextSZ(~,~)
        set(nextSZ,'Enable','off');
        drawnow;
        autoNext=1;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        tmFlag=0;set(f1,'Visible','off');
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
        Y_=LUT_previous(:,4);%
        Y_(indTopN(doneFlag==1))=LUTtopN_current(doneFlag==1,4);
        idx=find(ismember(Y_,'Seizure')==1);
        if isempty(idx)
            choice=questdlg('No seizure found!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
        else          
            tmp_=find((idx-jjj)>0,1,'first');
            if isempty(tmp_)
                choice=questdlg('No new seizure!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
            else
                jjj=idx(tmp_);
                uiresume(f);
            end
        end 
    end

    function fcn_prevSZ(~,~)
        set(prevSZ,'Enable','off');
        drawnow;
        tmFlag=0;
        autoNext=1;set(autoNextH,'value',autoNext);
        set(f1,'Visible','off');
        queryMode=0;
        set(queryModeH,'value',queryMode);
        Y_=LUT_previous(:,4);
        Y_(indTopN(doneFlag==1))=LUTtopN_current(doneFlag==1,4);
        idx=find(ismember(Y_,'Seizure')==1);
        if isempty(idx)
            choice=questdlg('No seizure found!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
        else
            tmp_=find((idx-jjj)<0,1,'last');
            if isempty(tmp_)
                choice=questdlg('No old seizure!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
            else
                jjj=idx(tmp_);
                uiresume(f);
            end
        end 
    end

    function fcn_tm(~,~)
        set(tm,'Enable','off');
        autoNext=0;set(autoNextH,'value',autoNext);
        queryMode=0;set(queryModeH,'value',queryMode);
        SeqMode='Sequential off';set(seqModePop,'value',1);
        ModMode='Model off';set(modModePop,'value',1);
        
        tmp=get(bg,'SelectedObject');p=tmp.String;          
        if strcmp(p,'None')
            tmFlag=-1;
            choice=questdlg('Please pick a sample and label it 1st!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    set(0,'CurrentFigure',f)
                    set(f1,'Visible','off');
                    uiresume(f);
            end            
        else
            tmFlag=1;            
            ip=regexp(p,']');
            y_tm=p(ip+2:end);        
            
            set(f1,'Visible','on');
            set(goBack,'Enable','on');
            for id_r=1:nP
                set(r1{id_r},'Enable','on')
            end
            set(ro1,'Enable','on');
            set(bg1,'SelectedObject',ro1)
            for iAx=1:length(Ax_spec1)
                set(Ax_spec1{iAx},'Visible','on');
            end
            
            set(0,'CurrentFigure',f1)
            indSeen_=unique([indSeen_;indTopN(doneFlag==1)]);
            indUnseen_=setdiff((1:nSamples)',indSeen_);     
            if isempty(indUnseen_)
                tmFlag=-1;
                choice=questdlg('All CP center got labeled!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        set(0,'CurrentFigure',f)
                        set(f1,'Visible','off');
                        uiresume(f);
                end                
            else
                bow_vec_unseen=bow_vec(indUnseen_,:);bow_vec_current=bow_vec(jjj,:);
                D_=fcn_distChiSq(bow_vec_unseen,bow_vec_current);
                [~,idx_]=min(D_);
                jjj1=indUnseen_(idx_);
                iEpoch1=LUT_previous{jjj1,1};
                ii1=LUT_previous{jjj1,2};         
                
                tmp=load([dataDir,iEpoch1,'_',num2str(ii1),'.mat']);
                seg1=tmp.SEG;sdata1=tmp.Sdata;sfreqs1=tmp.sfreqs;tc1=2*ii1-1;timeStampo1=datetime('2001-01-01 00:00:00');
                eeg1=seg1(1:19,:);eeg1(eeg1>300)=300;eeg1(eeg1<-300)=-300;
                ekg1=seg1(20,:);
                w=14;
                tto1=tc1*Fs-(w/2)*Fs+1;tt11=tc1*Fs+(w/2)*Fs;tt_1=tto1:tt11;
                gap=NaN(1,size(eeg1,2));
                switch montage
                    case 'L-Bipolar'
                        seg=fcn_bipolar(eeg1);
                        seg_disp1=[seg(1:4,:);gap;seg(5:8,:);gap;seg(9:12,:);gap;seg(13:16,:);gap;seg(17:18,:);gap;ekg1];
                        channel_withspace1=channel_withspace_bipolar;
                    case 'Average' 
                        seg=eeg-repmat(mean(eeg,1),size(eeg,1),1);
                        seg_disp1=[seg(1:8,:);gap;seg(9:11,:);gap;seg(12:19,:);gap;ekg1];
                        channel_withspace1=channel_withspace_average;
                    case 'Monopolar'
                        seg=eeg;
                        seg_disp1=[seg(1:8,:);gap;seg(9:11,:);gap;seg(12:19,:);gap;ekg1];
                        channel_withspace1=channel_withspace_monopolar;
                end
                M1=size(seg_disp1,1);
                DCoff1=repmat(flipud((1:M1)'),1,size(seg_disp1,2));
                timeStamps1=datestr(timeStampo1+seconds(round(tto1/Fs:2:(tt11+1)/Fs)),'hh:MM:ss');
                ekg_1=seg_disp1(end,:);ekg_1=(ekg_1-nanmean(ekg_1))/(eps+nanstd(ekg_1));     
                set(f1,'CurrentAxes',Ax_EEG1);cla(Ax_EEG1)
                hold(Ax_EEG1,'on')
                    for iSec=1:round((tt11-tto1+1)/Fs)
                        ta=tto1+Fs*(iSec-1);
                        line(Ax_EEG1,[ta ta],[0 M1+1],'linestyle','--','color',[.5 .5 .5])
                    end
                    plot(Ax_EEG1,tt_1,zScale*seg_disp1(1:end-1,:)+DCoff1(1:end-1,:),'k');         
                    plot(Ax_EEG1,tt_1,-.2*ekg_1+DCoff1(end,:),'m');
                    set(Ax_EEG1,'ytick',1:M1,'yticklabel',fliplr(channel_withspace1),'box','on','ylim',[0 M1+1],'xlim',[tto1 tt11+1],'xtick',round(tt_1(1):2*Fs:(tt_1(end)+Fs)),'xticklabel',timeStamps1)
                    
                    dt=tt11-tto1+1;a=round(dt*4/5);
                    xa1=tto1+[a a+Fs-1];xa2=tto1+[a a];ya1=[3 3];ya2=ya1+[0 100*zScale];
                    text(Ax_EEG1,xa1(1)-.7*a/10,mean(ya2),'100\muV','Color','r','FontSize',12);text(Ax_EEG1,mean(xa1)-0.3*a/10,2.5,'1 sec','Color','r','FontSize',12);
                    line(Ax_EEG1,xa1,ya1,'LineWidth',2,'Color','r');line(Ax_EEG1,xa2,ya2,'LineWidth',2,'Color','r');
                hold(Ax_EEG1,'off')
 
                winSize=10;
                to1=(tc1-1)-(winSize*60/2)+1;t11=(tc1-1)+(winSize*60/2);
                S_x1=to1:stepSize:t11;S_y1=sfreqs1;
                for iReg=1:4
                    set(f1,'CurrentAxes',Ax_spec1{iReg});cla(Ax_spec1{iReg})
                    spec=sdata1{iReg};
                    colormap jet;axis(Ax_spec1{iReg},'xy');           
                    hold(Ax_spec1{iReg},'on')
                        imagesc(Ax_spec1{iReg},S_x1,S_y1,pow2db(spec),col);
                        plot(Ax_spec1{iReg},[tc1  tc1],[S_y1(1) S_y1(end)],'k--','linewidth',1);
                    hold(Ax_spec1{iReg},'off')
                    
                    ylabel(Ax_spec1{iReg},'Freq (Hz)');
                    yticklabels=get(Ax_spec1{iReg},'yticklabel');yticklabels{end}=spatialRegs{iReg};
                    set(Ax_spec1{iReg},'yticklabel',yticklabels,'xtick',[],'box','on','ylim',[S_y(1) S_y(end)],'xlim',[tc1-winSize/2*60,tc1+winSize/2*60])      
                end
                
                tt_1=timeStampo1+seconds(to1:5*60:(t11+1));
                hold(Ax_spec1{4},'on')
                    set(Ax_spec1{4},'xtick',to1:5*60:(t11+1),'xticklabel',datestr(tt_1,'hh:MM:ss'),'box','on','xlim',[tc1-winSize/2*60,tc1+winSize/2*60]);
                hold(Ax_spec1{4},'off')  
                
                set(f1,'CurrentAxes',Ax_spec1{end});cla(Ax_spec1{end})
                plot(Ax_spec1{end},tc1,0,'rv','markersize',8,'MarkerFaceColor','r');
                set(Ax_spec1{end},'xtick',[],'xlim',get(Ax_spec1{end-1},'xlim'),'ylim',[0 .5]) 
                axis(Ax_spec1{end},'off');
            end
        end
        uiresume(f);
    end

    function fcn_done(~,~)
        tmFlag=0;
        autoNext=1;set(autoNextH,'value',autoNext);
        set(f1,'Visible','off');

        idxSeen=unique([indSeen_;indTopN(doneFlag==1)]); 
        Y_=LUT_previous(:,4);
        Y_(indTopN(doneFlag==1))=LUTtopN_current(doneFlag==1,4);
        LUT=LUT_previous(:,1:3);
        LUT(:,3)=Y_;
        
        timeStamp=datestr(now,'yyyy-mmm-dd HH:MM:ss.FFF');
        newLUTname=[targetDir,'ActiveLearning_LUT_final.mat'];
        idxReal_idx=cell2mat(idxReal(:,1));
        [~,idx]=unique(idxReal_idx,'last');
        idxReal=idxReal(idx,:);        
        
        save(newLUTname,'LUT','Y_model','Xvisual','bow_vec','idxSeen','timeStamp','rater','idxReal');
        
        choice=questdlg('All set!','Warning','Ok','Ok');
        switch choice
            case 'Ok'
                uiresume(f);
        end
    end

    function D=fcn_distChiSq(X,Y)
        mm=size(X,1);nn=size(Y,1);
        mOnes=ones(1,mm);D=zeros(mm,nn);
        for i=1:nn
            yi=Y(i,:);
            yiRep=yi(mOnes,:);
            s=yiRep+X;
            dd=yiRep-X;
            D(:,i)=sum(dd.^2./(s+eps),2);
        end
        D=D/2;
    end
end
