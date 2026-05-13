function [slope, t_out, freqs] = sm_wavelet_1overf(x, fs, step_sec, freq_range, n_freqs, nCycles)
% WAVELET_1OVERF  Memory-efficient wavelet-based 1/f slope estimation
%
% Inputs:
%   x          - time series vector (1 x N or N x 1)
%   fs         - sampling rate (Hz)
%   step_sec   - time resolution for slope (seconds)
%   freq_range - [fmin fmax] frequency range in Hz
%   n_freqs    - number of log-spaced frequencies for fitting
%   nCycles    - number of cycles per Morlet wavelet (time/freq tradeoff)
%
% Outputs:
%   slope      - 1/f slope over time
%   t_out      - time vector corresponding to slope estimates
%   freqs      - frequencies used for slope fitting

%% Ensure column vector
x = x(:);
N = length(x);

%% Log-spaced frequencies
freqs = logspace(log10(freq_range(1)), log10(freq_range(2)), n_freqs);
logf = log10(freqs);

%% Time points for slope computation
step_samp = round(step_sec * fs);
time_idx = 1:step_samp:N;
t_out = (time_idx - 1) / fs;

slope = nan(1, length(time_idx));

%% Compute slope per time point
for ti_idx = 1:length(time_idx)
    ti = time_idx(ti_idx);
    power = nan(n_freqs,1);

    for fi = 1:n_freqs
        f = freqs(fi);
        sigma_t = nCycles / (2*pi*f);
        t_wav = -3*sigma_t : 1/fs : 3*sigma_t;

        idx = ti + round(t_wav*fs);
        valid_idx = idx > 0 & idx <= N;

        x_seg = x(idx(valid_idx));
        wave_seg = exp(2*1i*pi*f*t_wav(valid_idx)) .* exp(-t_wav(valid_idx).^2 / (2*sigma_t^2));

        x_seg = x_seg(:);
        wave_seg = wave_seg(:);

        conv_val = sum(x_seg .* conj(wave_seg));
        power(fi) = abs(conv_val)^2;   % wavelet power
    end

    % Fit slope in log-log space
    valid = power > 0;
    if sum(valid) >= 5
        p = polyfit(logf(valid), log10(power(valid)), 1);
        slope(ti_idx) = p(1);
    end
end
end
