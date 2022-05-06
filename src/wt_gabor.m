
% load raw signal
matfile = load("/Users/steve/Documents/code/unm/IHKApy/src/raw_signal.mat");
signal = matfile.signal;
fs=16.0; % arbitrary
freqlist=[0.5 1.0 2.0 4.0 8.0]; % log spacing
wt = wt_matlab(signal,fs,freqlist,"Gabor",5);

% save the transformed array
save("wt_ml.mat","wt")


% wt_matlab aka awt_freqlist
% SF: this function was previously in the buzcode repository in /externalPackages/ directory
function [wt,freqlist,psi_array] = wt_matlab(x,Fs,freqlist,type,xi)
%   awt_freqlist  analytical wavelet transform, where one can specify the list of desired frequencies 
%   
%   [wt,freqlist,psi_array] = awt_freqlist(x,Fs,freqlist,type,xi)
%
%   Inputs:
%       x           the signal to be analyzed
%       Fs          the sampling frequency
%       freqlist    list of frequencies at which to compute wt (optional)
%                   (or set to [] for automatic definition)
%       type        type of wavelet to use (Gabor, Lusin or Sombrero)
%       xi          the number of oscillations parameter
%   Outputs:
%       wt: time-frequency image
%       freqlist: useful in case of automatic definition of frequencies
%       psi_array : array of analysis functions (complex values)
%
%  Maureen Clerc, Christian Benar, october 2007
%  modified from awt from wavelab toolbox

% History of changes
% 1/11/2007: (cgb) psi_array: output in complex 
% 3/06/2008: (cgb) init of psi_array to size of wt

n = length(x);
sigma2 = 1;
x = x(:);
omega = [(0:n/2) (-ceil(n/2)+1:-1)].*Fs/n; % CGB added ceil: TO BE CHECKED
omega = omega(:);

fftx = fft(x);

% to compute min and max center frequency values:
%tolerance = 1.5; % such that an acceptable "zero" is exp(-tolerance^2/2)
tolerance = 0.5; % cgb
%tolerance = 1; % cgb

if nargin<2 
    Fs = 1;
end
if nargin<4
    type = 'Gabor';
end
if nargin<5
    xi = 5; % only useful for Gabor
end

% compute min and max valid center frequencies, according to wavelet type
mincenterfreq = 2*tolerance*sqrt(sigma2)*Fs*xi./n;
maxcenterfreq = Fs*xi/(xi+tolerance/sqrt(sigma2));

if  (nargin<3 || isempty(freqlist) )
    nvoice = 12;
    freqlist= 2.^(log2(mincenterfreq):1/nvoice:log2(maxcenterfreq));
end

s_array = xi./freqlist;
minscale = xi./maxcenterfreq;
maxscale = xi./mincenterfreq;

nscale = length(freqlist);
wt = zeros(n,nscale);
scaleindices=find(s_array(:)'>=minscale & s_array(:)'<=maxscale);
%psi_array=zeros(n,length(scaleindices));
psi_array=zeros(n,nscale);
for kscale=scaleindices
    s=s_array(kscale);

    freq =  (s .* omega  - xi);
    Psi = realpow(4.*pi.*sigma2,1/4)*sqrt(s) *exp(-sigma2/2*freq.*freq);
            

    wt(1:n,kscale) = ifft(fftx.*Psi);
    psi_array(:,kscale)=ifft(Psi);
end

end

