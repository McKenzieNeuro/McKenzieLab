function [EEG_clean, isBadChannels] = fcn_cleaningPipeline(EEG, EKG, verbose)
    
    tmp = load('LOC_18channels_polar.mat');
    Vxy = cell2mat(tmp.LUT(:,2)); 

    tmp = load('BC_LUT_v2.mat');
    BC_LUT = tmp.LUT;
    K = 10; % maintain 99% train senstivity %

    % EEG: mono-polar 19 channels %
    Fs = 200;
    win = 10*Fs; % 10sec %
   
    [M, N] = size(EEG);
    nWin = ceil(N/win);
    
    EEG_tmp = zeros(M, nWin*win);
    EKG_tmp = zeros(1, nWin*win);
    
    EEG_tmp(:, 1:N) = EEG; 
    EKG_tmp(:, 1:N) = EKG;
    
    EEG_clean = NaN(M-1, nWin*win);
    isBadChannels = NaN(nWin, M-1);
    for i = 1:nWin
        a = (i-1)*win+1;
        b = i*win;
        
        eeg = EEG_tmp(:, a:b);
        ekg = EKG_tmp(:, a:b);
        
        
        if isnan(mean(ekg)) || mean(abs(ekg)) == 0 % NaNs or all 0s %
            ekg = mean(eeg, 1);
        end
        eeg = fcn_Bipolar(eeg);
        
        % step2: detect bad channels %
        isBad = fcn_isBadChannel_v2(eeg, ekg, Fs, BC_LUT, K);

        % step3: handle bad channel by 2D interpolation %
        eeg_clean = fcn_fixBadChannel(eeg, isBad, Vxy);

        % step4: export %
        EEG_clean(:, a:b) = eeg_clean;
        isBadChannels(i, :) = isBad;
        
        
        % step5: plot for sanity check - bad in red, otherwise black %
        if verbose
            figure('units','normalized','outerposition',[0.5 0.5 .5 .5]);
            hold on
            for j = 1:18
                if isBad(j)
                    plot((1:size(eeg, 2))/Fs, eeg(j,:)/200+18-j, 'color', [.7 .7 .7]);
                    plot((1:size(eeg, 2))/Fs, eeg_clean(j,:)/200+18-j, 'r');
                else
                    plot((1:size(eeg, 2))/Fs, eeg(j,:)/200+18-j, 'k');
                end

            end
            box on
            ylim([-1 18])
            hold off
            close
        end
    end
    
    % adjust in length %
    EEG_clean = EEG_clean(:, 1:N);

end