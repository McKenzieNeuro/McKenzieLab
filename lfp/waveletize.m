
function [POWER,PC,frex] = waveletize(DataIn,srate,numFreq,lowCut,hiCut,space,Par,varargin)


if ~isempty(varargin)
    
    freqBand = varargin{1};
else
    freqBand = 1:numFreq;
end


    % Setup Wavelet Params
    if space == 0
        frex = logspace(lowCut,log10(hiCut),numFreq);
    else
        frex = linspace(lowCut,hiCut,numFreq);
    end
    s = logspace(log10(3),log10(10),numFreq)./(2*pi*frex);
    t = -2:1/srate:2;
    WAVELET_TIME = repmat([1,size(DataIn,1)],1);

    % Do time-freq

    if size(DataIn,2) > 1
        DataIn = bsxfun(@minus,DataIn,mean(DataIn,1));
    end

    % Definte Convolution Parameters
    dims = size(DataIn);
    n_wavelet = length(t);
    n_data = dims(1)*dims(2);
    n_convolution = n_wavelet+n_data-1;
    n_conv_pow2 = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;

    % get FFT of full range of data
    EEG_fft = fft(reshape(DataIn,1,n_data),n_conv_pow2); 

    POWER = zeros(length(freqBand),WAVELET_TIME(1,2));
    PC = zeros(length(freqBand),WAVELET_TIME(1,2));
    
    if Par == 0
        for iFig = freqBand

            wavelet = fft(exp(2*1i*pi*frex(iFig).*t).*exp(-t.^2./(2*(s(iFig)^2))),n_conv_pow2 );
            % convolution
            EEG_conv = ifft(wavelet.*EEG_fft);
            EEG_conv = EEG_conv(1:n_convolution);
            EEG_conv = EEG_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
            EEG_conv = reshape(EEG_conv,dims(1),dims(2));
            % Get power
            POWER(iFig,:) = abs(mean(EEG_conv,2)).^2;

            % Get ITPC by condi (different time windows)
            PC(iFig,:) = abs(mean(exp(1i*(angle(EEG_conv(:,:)))),2));
        end
    else
        parfor iFig = freqBand

            wavelet = fft(exp(2*1i*pi*frex(iFig).*t).*exp(-t.^2./(2*(s(iFig)^2))),n_conv_pow2 );
            % convolution
            EEG_conv = ifft(wavelet.*EEG_fft);
            EEG_conv = EEG_conv(1:n_convolution);
            EEG_conv = EEG_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
            EEG_conv = reshape(EEG_conv,dims(1),dims(2));
            % Get power
            POWER(iFig,:) = abs(mean(EEG_conv,2)).^2;

            % Get ITPC by condi (different time windows)
            PC(iFig,:) = abs(mean(exp(1i*(angle(EEG_conv(:,:)))),2));
        end
    end
end
        