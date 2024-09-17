function [Qa] = InterpolationUniformDt(Qt,Qp,Klist,ta)

    % Takes trials corresponding to time vectores with different values of
    % dt along the trial and across trials and converts them to a uniform
    % value of dt (both along and across trials)

    % Klist: list of trial numberes to convert
    
    
    Pa = zeros(1,length(ta));
    Qa = zeros(length(ta),length(Klist));
    for k=1:length(Klist)
        [Qtmax,jmax] = max(Qt(:,Klist(k)));
        tc = Qt(1:jmax,Klist(k));
        Pc = Qp(1:jmax,Klist(k));
        Pa(1) = Pc(1);
        for j=2:length(ta)        
            jcant = find(tc<ta(j),1,'last');
            tcant = tc(jcant);
            jcpost = find(tc>ta(j),1,'first');
            tcpost = tc(jcpost);
            Pcant = Pc(jcant);
            Pcpost = Pc(jcpost);
            Pa(j) = Pcant+(ta(j)-tcant)/(tcpost-tcant)*(Pcpost-Pcant);
        end
        Qa(:,k) = Pa;   
        k 
    end