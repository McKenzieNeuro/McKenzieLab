
function Y_f = sm_notch_filter(y,fs,notch)
          
fn = fs/2;              % Nyquist frequency
freqRatio = notch/fn;      % ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       % width of the notch

% Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
notchPoles = (1-notchWidth) * notchZeros;

%figure;
%zplane(notchZeros.', notchPoles.');

b = poly( notchZeros ); %  Get moving average filter coefficients
a = poly( notchPoles ); %  Get autoregressive filter coefficients

%figure;
%freqz(b,a,32000,fs)

% filter signal x
Y_f = filter(b,a,y);
end