clear all
cd('R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220630-112150') % good NO
cd('R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\LinearTrack\NE2h2-220701-090608') %con
%%
cd('R:\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220630-125426')
cd('R:\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220720-090530') % pretty good
cd('R:\DANEHippocampalResponse\NE2h6\Linear Track\NE2h6-220722-091944') % control
%cd('R:\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220721-082011')
%cd('R:\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220720-093937')
%cd('R:\DANEHippocampalResponse\NE2h9\Linear Track\NE2h9-220720-102405')
fils = getAllExtFiles(pwd,'avi',0);%faster to read from local
vid = VideoReader(fils{1});
%vid = vision.VideoFileReader(fils{1});
load('sessiondata.mat')

fs = sessiondata.neural.fs_neural;
tim_in = sessiondata.contextEntry{4,2};

[ix,early,late,ts] = sm_getIndicesAroundEvent(tim_in,10,125,fs,length(sessiondata.neural.signal_DFoF));
k = gaussian2Dfilter([ 1 fs*10],fs);
sig = nanconvn(sessiondata.neural.signal_DFoF,k);



ts = ts+tim_in;
ts_all = (1:length(sig))/fs;

nix = round(10*fs);
%%
close all
v = VideoWriter("newfile5.avi");
open(v)



nix = round(10*fs);

for i = 1:50:length(ts)
    tic
    t= ts_all(ix(i)-nix:ix(i)+nix);
    [~,b] = bestmatch(ts(i),sessiondata.behavior.ts_video);
    
    
  im = image('XData',[1 20347],'CData',read(vid,b));
hold on
rectangle('Position',[0 -200 20347 300],'FaceColor','w')
plot(85*sig(ix(i)-nix:ix(i)+nix)-75,'k')
plot([1 20347],[0 0]-75,'--','color','k')
plot([20347 20347]/2,[-200 100],'--','color','k')

plot([0 (20347/20)]+500,[0 0],'k')
plot([0 0]+500,[0 85],'k')
set(gca,'ydir','normal')
hold off
 
  ylim([-200 460])
  xlim([0 20347])
  axis off
 frame = getframe(gcf); %get frame
 writeVideo(v,frame)

toc
  clf

end
close(v)
