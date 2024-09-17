function [IaF,NegValuesPos,xm,dm]=InstFreq_v2(x,dt,typeEnv,fig,sm,verbose)
%
% Instantaneous angular frequency IaF of the signal x.
%
% Output
%
%    IaF = Instantaneous angular frequency of the given signal x
%    NegValuesPos = position of the negative values
%    xm = envelope of the signal x
%    dm = envelope of the derivative of the signal
%
% Input
%
%   x     = signal
%
%   dt    = time interval
%
% typeEnv = 'linear' the envelope is a piecewise linear function
%                     connecting abs of min and max of f
%            ('cubic')  the envelope is a piecewise cubic function
%                     connecting abs of min and max of f
%            'PERFcubic' the envelope is a piecewise cubic function connecting
%                    abs of min and max of f and everywhere |f|<=Env.
%            'EnoCubic' Essentially Non-Oscillatory (ENO) technique
%
%   fig    = (0) no figures displayed
%             1 figures are generated
%
%   sm     = (0) no smoothing of the Instantaneous angular frequency
%             1  smoothing
%
%   verbose = (0) no output displayed,
%              1  otherwise
%
%   See also SMOOTH, CONSENV, DETJUMP_V03.
%
% This software is an alpha version, please report any bug, comment or
% suggestion to cicone@math.gatech.edu
%
% References
%    A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for
%      Signal Decomposition and Instantaneous Frequency analysis'.
%      to appear in ACHA Special Issue �Sparse Representations with
%      Applications in Imaging Science, Data Analysis and Beyond�
%      Preprint ArXiv http://arxiv.org/abs/1411.6051
if nargin<3; typeEnv='cubic'; end
if nargin<4; fig=0; end
if nargin<5; sm=0; end
if nargin<6; verbose=0; end
if fig
    figure;
    plot(x,'k');
    title('Given signal')
end
xm=ConsEnv(x,typeEnv,fig);
notNull=xm~=0;
x(notNull) = x(notNull) ./xm(notNull) ;
if fig
    figure;
    plot(x,'k');
    title('Scaled signal')
end
d=(x(3:end)-x(1:end-2))./(2*dt);
if fig
    figure;
    plot(d,'k');
    title('Signal derivative')
end
N=length(d);
dm=ConsEnv(d,typeEnv,fig);
notNull=dm~=0;
d(notNull) = d(notNull)./dm(notNull);
if fig
    figure;
    plot(dm/2/pi,'k');
    title('Derivative envelope')
end
if fig
    figure;
    plot(d,'k');
    title('Scaled derivative')    
end
xx=x(2:end-1);
if fig
    figure
    plot(xx,d)
    title('Rotation after scaling')
end
phase=atan2(xx,d);
if fig
    figure
    plot(phase);
    title('Original phase')
end
countNeg=0;
NegValuesPos=[];
for i=2:N-1
    if phase(i-1)>phase(i)
        if phase(i-1)-phase(i)>pi
            phase(i:N)=phase(i:N)+2*pi;
        else
            phase(i)=phase(i-1);
            countNeg = countNeg +1;
            NegValuesPos=[NegValuesPos i];
        end
    end
end
if countNeg>0 && verbose>0
    fprintf('         WARNING: Negative phase in %1.0d points \n',countNeg)
end
if fig
    figure
    plot(phase);
    title('Phase after shifting')
end
IaF = (phase(3:end)-phase(1:end-2))./(dt*4*pi);
if fig
    figure
    plot(IaF);
    title('Instantaneous angular frequency')
end
if sm
    IaF = smooth(IaF);
    
    if fig
        figure
        plot(IaF);
        hold all;
        title('Smooth instantaneous angular frequency')
    end
end
end
function[T1] = smooth(T)
%Take mean of envelopes
x = T;
N=length(x);
for i=1:N
    if isnan(x(i))==1
        x(i) =( x(i+1)+x(i-1) )/2;
    end
end
d=diff(x);
maxmin = []; % to store the optima (min and max without distinction so far)
for i=1:N-2
    if d(i)==0                        % we are on a zero
        maxmin = [maxmin, i];
    elseif sign(d(i))~=sign(d(i+1))   % we are straddling a zero so
        maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
    end
end
% k=length(maxmin);
% divide maxmin into maxes and mins
if maxmin(1)>maxmin(2)              % first one is a max not a min
    maxes = maxmin(1:2:length(maxmin));
    mins  = maxmin(2:2:length(maxmin));
else                                % is the other way around
    maxes = maxmin(2:2:length(maxmin));
    mins  = maxmin(1:2:length(maxmin));
end
k1=length(maxes);   k2=length(mins);
XX1=[1,maxes,N];
YY1=[x(maxes(1)),x(maxes),x(maxes(k1))];
XX2=[1,mins,N];
YY2=[x(mins(1)),x(mins),x(mins(k2))];
maxenv=interp1(XX1,YY1,1:N,'cubic');
minenv = interp1(XX2, YY2,1:N,'cubic');
T1 = (maxenv+minenv)/2;
end
function [Env,maxmin]=ConsEnv(f,typeEnv,fig)
%
% Construct the envelope Env of the signal f.
%
% Output
%
%    Env = envelope
%    maxmin = position of Max and Min of the given signal f
%
% Input
%
%   f      = signal
% typeEnv = ('linear') the envelope is a piecewise linear function
%                     connecting abs of min and max of f
%            'cubic'  the envelope is a piecewise cubic function
%                     connecting abs of min and max of f
%            'PERFcubic' the envelope is a piecewise cubic function connecting
%                    abs of min and max of f and everywhere |f|<=Env.
%            'ENOcubic' Essentially Non-Oscillatory (ENO) technique
%
%   fig    = (0) no figures displayed
%             1 figures are generated
%
if nargin<2; typeEnv='linear'; end
if nargin<3; fig=0; end
epsm=10^-14;
N=length(f);
d=diff(f);
maxmin = []; % to store the optima (min and max without distinction)
for i=1:N-2
    if sign(d(i))~=sign(d(i+1))
        maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
    end
end
if strcmp(typeEnv,'PERFcubic')
    
    %% PERFcubic interpolation
    
    Env=interp1([1 maxmin N],[max(abs(f(maxmin(1))),abs(f(1))) abs(f(maxmin)) max(abs(f(maxmin(end))),abs(f(end)))],1:N,'cubic');
    
    if fig
        figure
        plot(f,'k')
        axis([1 N 1.1*min(f) 1.1*max(f)]) %-3*10^-3 3*10^-3])
        hold on
        plot(Env,'r')
        plot(-Env,'g')
        plot(zeros(1,length(f)),'b')
        title('Initial Envelope')
        legend('Signal','Env','-Env')
    end
    
    loop1=0;
    k=0;
    if not(isempty(find(abs(f)-Env>epsm)))
        loop1=1;
    end
    
    while loop1 && k<50
        loop1=0;
        k=k+1;
        fprintf('Step k = %1.0d\n',k)
        
        % adding back missing almost-max/min
        fs=max(abs(f)-Env,epsm);
        ds=diff(fs);
        imaxmin = [];
        
        for i=1:N-2
            if sign(ds(i))~=sign(ds(i+1)) && sign(ds(i))>sign(ds(i+1))
                imaxmin = [imaxmin, i+1];
            end
        end
        
        fprintf(' %1.0d new points added\n\n',length(imaxmin))
        maxmin=sort([maxmin imaxmin]);
        
        Env=interp1([1 maxmin N],[max(abs(f(maxmin(1))),abs(f(1))) abs(f(maxmin)) max(abs(f(maxmin(end))),abs(f(end)))],1:N,'cubic');
        
        if not(isempty(find(abs(f)-Env>epsm)))
            loop1=1;
        end
    end
    
    if not(isempty(find(abs(f)-Env>epsm)))
        disp('Doesn''t work!')
        disp('Try increasing k')
        return
    end
    
elseif strcmp(typeEnv,'ENOcubic')
    
    XX = [1, maxmin,N];
    YY = [max(abs(f(maxmin(1))),abs(f(1))),abs(f(maxmin)),max(abs(f(maxmin(end))),abs(f(end)))];
    
    i=1;
    XXX1=XX(1);
    YYY1=YY(1);
    
    Env=zeros(1,N);
    
    while i< length(YY)
        
        % no jump detected
        if abs( YY(i) - YY(i+1) ) < 0.4* max(YY(i), YY(i+1)) %threshold 0.4
            XXX1=[XXX1, XX(i+1)];
            YYY1=[YYY1, YY(i+1)];
        else % jump detected
            % this code can handle only one jump between two extreme points
            % of a given signal
            index = DetJump_v03(f,XX(i),XX(i+1));
            if index==XX(i) % we do no add any intermidiate point, observe that index cannot be ever = XX(i+1)
                % XXX1=XXX1;
                % YYY1=[YYY1];
            else % we need to add an intermidiate point
                XXX1=[XXX1,index];
                YYY1=[YYY1, YY(i)];
            end
            Env(XXX1(1):XXX1(end)) = interp1(XXX1, YYY1, XXX1(1):index,'cubic');
            if index+1==XX(i+1) % we do no add any intermidiate point
                XXX1=XX(i+1);
                YYY1=YY(i+1);
            else
                XXX1=[index+1, XX(i+1)];
                YYY1=[YY(i+1), YY(i+1)];
            end
        end
        i=i+1;
    end
    Env(XXX1(1):XXX1(end)) = interp1(XXX1, YYY1, XXX1(1):XXX1(end),'cubic');
    
else  % linear or cubic
    XX = [1, maxmin,N];
    YY = [max(abs(f(maxmin(1))),abs(f(1))),abs(f(maxmin)),max(abs(f(maxmin(end))),abs(f(end)))];
    
    Env = interp1(XX, YY,1:N,typeEnv);
end
if fig
    figure
    plot(f,'k')
    axis([1 N 1.1*min(f) 1.1*max(f)]) %-3*10^-3 3*10^-3])
    hold on
    plot(Env,'r')
    plot(-Env,'g')
    plot(zeros(1,length(f)),'b')
    title('Envelopes')
    legend('Signal','Env','-Env')
end
end
function index = DetJump_v03(f,I1,I2)
%figure;plot(f)
d= diff(f(I1:I2));
[dvalue,dpos]=sort(abs(d));
%figure;plot(I1+0.5:I2-0.5,d)
dd= diff(d);
%figure;plot(I1+1:I2-1,dd)
[ddvalue,ddpos]=sort(abs(dd));
prompt = ['Default values for identifying a jump in the signal and its derivative are\n'...
    ' sJump = 0.25\n dJump = 0.25\n\n Do you want to change them? Y/N [N]: '];
str = input(prompt,'s');
if isempty(str)
    str = 'N';
end
if str=='N'
    sJump=0.25;
    dJump=0.25;
else
    sJump=input('Insert a value < 1 for sJump\n << ');
    dJump=input('Insert a value < 1 for dJump\n << '); 
end
if (ddvalue(end)-ddvalue(end-1))/ddvalue(end)<=dJump && (dvalue(end)-dvalue(end-1))/dvalue(end)<=sJump
    % the difference between the biggest jump and the almost biggest one is less than sJump times the biggest one
    disp('More than one jump in the signal and its derivative, or no jump at all!')
    disp('The current algorithm can handle one and only one jump per interval')
    index=I1;
elseif (ddvalue(end)-ddvalue(end-1))/ddvalue(end)>dJump && (dvalue(end)-dvalue(end-1))/dvalue(end)<=sJump
    % unique jump in dd
    if ddpos(end)==1
        index=ddpos(end)-1+I1;
    elseif ddpos(end)==length(dd)
        index=ddpos(end)+I1;
    else
        index=ddpos(end)+I1;
    end
elseif (ddvalue(end)-ddvalue(end-1))/ddvalue(end)<=dJump && (dvalue(end)-dvalue(end-1))/dvalue(end)>sJump
    % unique jump in d
    index=dpos(end)-1+I1;
else
    [value,pos]=max([(dvalue(end)-dvalue(end-1))/dvalue(end),(ddvalue(end)-ddvalue(end-1))/ddvalue(end)]);
    if pos==2
        if ddpos(end)==1
            index=ddpos(end)-1+I1;
        elseif ddpos(end)==length(dd)
            index=ddpos(end)+I1;
        else
            index=ddpos(end)+I1;
        end
    else
        index=dpos(end)+I1-1;
    end
end
end
