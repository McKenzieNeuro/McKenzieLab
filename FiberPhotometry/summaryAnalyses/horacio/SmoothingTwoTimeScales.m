function [Qasmoothcomb] = SmoothingTwoTimeScales(Qa,span1,span2,tbdry,k,dt)

    % Takes the signal Qa
    % Smooth it using span1 (short time scale)
    % Computes the minimum between t=0 and t=tbdry
    % Smooth the smoothed signal between tbdry and tend using span2 (long
    % time scale)
    % Computes the combined smoothed signal
    
    Qasmooth = smooth(Qa(:,k),span1);
    [Qauxmin,jmin] = min(Qasmooth(1:floor(tbdry/dt)));
    Qasmooth1 = smooth(Qasmooth(jmin:end),span2);
    Qasmoothcomb = [Qasmooth(1:jmin-1)' Qasmooth1'];