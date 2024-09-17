function SEG = fcn_parseData(dataSize, seg, Fs, t_left, t_right, win_time)

    M = dataSize(1);
    N = dataSize(2);
    
    
   % ii: index of 2sec % 
   
    %t_center = (2*(ii-1)+1)*Fs+1;           % in pts %
    %t_left  = t_center - (Fs*win_time/2);
    %t_right = t_center + (Fs*win_time/2)-1;

    % EEG 50 sec with this in the center %
    SEG = zeros(M, win_time*Fs);
        
    % reasoning %
    if t_left<1  && t_right<N
        disp(' Early')
        a = 1-t_left+1;
        b = size(SEG, 2);
        SEG(:, a:b) = seg;

    elseif t_left<1  && t_right>=N
        disp(' Short')
        a = 1-t_left+1;
        b = a+N-1;
        SEG(:, a:b) = seg;

    elseif t_left>=1 && t_right<N
        disp(' Regular')
        SEG = seg;

    elseif t_left>=1 && t_right>=N
        disp(' Late')
        a = 1;
        b = N-t_left+1;
        SEG(:, a:b) = seg;
    else
        keyboard
    end
end