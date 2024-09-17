function ses = GetTransAcrossTime(times,groups,varargin)


duration = .1;
binSize = 0.0008;
Fs = 1/30000;
win = 200;
win_inc = 5;
win_std = .015;
int_win = [0.0008 .0032];

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'binsize'
            binSize = varargin{i+1};
            
       
        case 'duration'
            duration = varargin{i+1};
       
            
        case 'win_std'
            win_std = varargin{i+1};
         
            
            
            
        case 'win'
            win = varargin{i+1};
           
            
            
        case 'win_inc'
            win_inc = varargin{i+1};
       
            
        case 'Fs'
            Fs = varargin{i+1};
         
        case 'int_win'
            int_win = varargin{i+1};
            
            
            
    end
end

    conv_w = round(win_std/binSize);
    intwini = round(((int_win(1):binSize:int_win(2))/binSize) + duration/binSize/2+1);
    
    [ses.ccg,ses.cnt,ses.t] = CCGTime(times,groups,'Fs',Fs,'win_inc',win_inc,'win',win,'duration',...
        duration,'binsize',binSize);
    ses.trans = nan(size(ses.ccg,2),size(ses.ccg,3),size(ses.ccg,4));
    ses.prob = nan(size(ses.ccg));
    for k = 1:size(ses.ccg,4)
        for n = 1:size(ses.ccg,2)
            for m = 1 :size(ses.ccg,3)
                if n~=m
                    cch = ses.ccg(:,n,m,k);
                    
                    [pvals,pred,qvals]=cch_conv(cch,conv_w);
                    ses.prob(:,n,m,k) = (cch-pred)./ses.cnt(:,n,m,k);
                    
                    ses.trans(n,m,k) = nansum(ses.prob(intwini,n,m,k));
                end
            end
        end
    end
    


ses.win = win;
ses.intwini  =intwini;
ses.int_win = int_win;
ses.duration = duration;
ses.binSize = binSize;
ses.win = win;
ses.win_inc = win_inc;
ses.conv_w = conv_w;
ses.win_std = win_std;

