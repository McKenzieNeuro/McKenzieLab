function SEG = fcn_parseData(dd, seg, Fs, t1, t2, ww)

    % Initialization  
    mm = dd(1); nn = dd(2);
    SEG = zeros(mm, ww*Fs);
        
    % Reasoning  
    if t1<1&&t2<nn
        a = 1-t1+1; b = size(SEG, 2);
        SEG(:, a:b) = seg;
    elseif t1<1&&t2>=nn
        a = 1-t1+1; b = a+nn-1;
        SEG(:, a:b) = seg;
    elseif t1>=1&&t2<nn
        SEG = seg;
    elseif t1>=1&&t2>=nn
        a = 1; b = nn-t1+1;
        SEG(:, a:b) = seg; 
    else
        keyboard       
    end
end
