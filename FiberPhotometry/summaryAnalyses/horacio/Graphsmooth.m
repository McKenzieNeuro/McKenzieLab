function Graphsmooth(Qt,Qp,k,span)

    [Qtmax,jmax] = max(Qt(:,k));
    tc = Qt(1:jmax,k);
    Pc = Qp(1:jmax,k);
    figure
    hold on
    plot(tc,Pc,'-b','linewidth',2)
    set(gca,'fontsize',24);
    xlabel('t [sec]');
    ylabel('S');
    
    Pcs = smooth(Pc,span);
    plot(tc,Pcs,'-r','linewidth',2);
