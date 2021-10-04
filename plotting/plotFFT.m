function [f,amp,ang]= plotFFT(y,Fs,pl)


L= length(y);
T = 1/Fs;                     % Sample time
% Length of signal
t = (0:L-1)*T;                % Time vector


NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
amp = 2*abs(Y(1:NFFT/2+1));
ang = angle(Y(1:NFFT/2+1));
if pl
    
    % Plot single-sided amplitude spectrum.
    plot(f,amp)
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end
end