function [data2, channels2, ch_miss] = fcn_preprocess(data1, channels1, Fs1, Fs2, lut_map, fc)

    % map channels to what we know of %
    idx = find(ismember(channels1, lut_map(:, 1)));
    
    data1b = double(data1(idx, :));
    channels1b = channels1(idx);   
    channels1c = channels1b;
    for i = 1:size(channels1b, 1)
        idx = find(ismember(lut_map(:, 1), channels1b{i}));
        channels1c{i} = lut_map{idx, 2};
    end

    % re-arrange %
    channels2 = lower({  'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';...
                         'Fz'; 'Cz';'Pz';...
                         'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'});        
    ch_miss = []; ch_keep = [];
    for i = 1:19
        ch = channels2{i};
        
        idx = find(ismember(channels1c, ch));
        if isempty(idx)
            idx = find(ismember(channels1c, [ch, '_alt']));
        end
        
        if isempty(idx)
            ch_miss = [ch_miss, i];
        else
            if length(idx)>1
                keyboard
            else
                ch_keep = [ch_keep; [i, idx(1)]];
            end
        end
    end
    
    data2 = zeros(19, size(data1b, 2));
    data2(ch_keep(:, 1), :) = data1b(ch_keep(:, 2), :);
    ch_miss = channels2(ch_miss);
    
    % EKGs %
    idx = find(cellfun(@isempty, regexpi(channels1c, 'EKG'))==0);
    if isempty(idx) % No EKG - zeros
        ekg = zeros(1, size(data1b, 2));
    else
        [ch, ii] = sort(channels1c(idx));
        ekg = data1b(idx(ii), :);
        
        idx = find(ismember(ch, 'ekg')==1);
        idx1 = find(ismember(ch, 'ekg1')==1);
        idx2 = find(ismember(ch, 'ekg2')==1);
        
        if ~isempty(idx) % EKG found
            ekg = ekg(idx(1), :);
            
        else
            if ~isempty(idx1) && ~isempty(idx2) % Both L&R
                ekg = ekg(idx1(1), :) - ekg(idx2(1), :);  
            else 
                if ~isempty(idx1) && isempty(idx2) % Only L
                    ekg = ekg(idx1(1), :);
                else % Only R
                    ekg = ekg(idx2(1), :);
                end
            end
        end
        
    end
    
    % resample
    Fs1 = round(Fs1);
    Fs2 = round(Fs2);  
    if isempty(data2)
        data2 = []; 
    else
        if Fs1~=Fs2
            data2 = resample([data2;ekg]', Fs2, Fs1)';
        else
            data2 = [data2;ekg];
        end
        channels2 = [channels2; 'EKG'];
    end
    
    % denoise
    [B1, A1] = butter(3, [fc-2.5, fc+2.5]/(Fs2/2), 'stop');
    [B2, A2] = butter(3, [.5, 40]/(Fs2/2));
    
    data2 = filtfilt(B1, A1, data2')';
    data2 = filtfilt(B2, A2, data2')';
    
end