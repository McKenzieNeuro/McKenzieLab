function out = convLPFilt(data,cutOff,sampleFreq,filtLength,downSampleRatio)

    % data = input data, processed column-wise
    % cutOff = lowpass filter cutoff frequency, must not be > sampleFreq/2
    % sampleFreq = sampling frequency of input data
    % filtLength = filter length, [] uses default
    % downSampleRatio = ratio to downsample filtered data, [] to not use

    cutOff = cutOff/(sampleFreq/2);

    if isempty(filtLength)
        filtLength = 1250; % is this the low pass?
    end

    % create filter kernel    
    n = -floor(filtLength/2):floor(filtLength/2); % create time base
    B = sinc(cutOff*n).*cutOff; % make kernel     

    for i = 1:size(data,2)

        P = numel(data(:,i));
        Q = numel(B);
        filtLength = P+Q-1;
        K = 2^nextpow2(filtLength);
    
        % do the convolution
        afft = fft(data(:,i),K);
        bfft = fft(B,K);
        c = ifft(afft(:).*bfft(:));
    
        range = [floor(length(B)/2)+1,filtLength-ceil(length(B)/2)+1];
        data(:,i) = c(range(1):range(2));
    end

    if ~isempty(downSampleRatio)
        data = downsample(data,downSampleRatio);
    end

    out = data;
end