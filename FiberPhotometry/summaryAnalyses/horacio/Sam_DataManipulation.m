clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];


KSE = 2;
    % 1: horacio.mat: Data set 2023-08-01
    % 2: horacio1.mat: Data set 2023-08-02

if KSE == 1
    load('horacio.mat')
    t = X;
    P = Y;


    M = 0;
    for j=1:length(t)-1
        if t(j) > t(j+1)        
            M = [M j];
        end    
    end

    Ntr = length(M)-1;
    Msze = diff(M);

    Qt = zeros(max(Msze),Ntr);
    Qp = zeros(max(Msze),Ntr);
    for k=1:Ntr
        Qt(1:Msze(k),k) = t(M(k)+1:M(k+1));
        Qp(1:Msze(k),k) = P(M(k)+1:M(k+1));
    end

    figure
    imagesc(Qt')
    set(gca, 'Ydir','normal')
    colorbar

    figure
    imagesc(Qp')
    set(gca, 'Ydir','normal')
    colorbar

    % Trials starting for t(1)>0 (exclude)
    Aux1 = find(Qt(1,:)>0);
    % Trials ending to soon or too short (exclude)
    Aux4 = find(Qt(11000,:)==0);
    Aux1 = union(Aux1,Aux4);
    % Trials having abnormal traces (sudden large peaks)
    Aux2 = ceil(find(P>10)/max(Msze));
    Aux3 = zeros(1);
    Aux3(1) = Aux2(1);
    cntaux3 = 1;
    for j=2,length(Aux2);
        if Aux2(j)>Aux3(cntaux3)
            cntaux3 = cntaux3+1;
            Aux3(cntaux3) = Aux2(j);
        end
    end
    % Combined abnormal trials
    Aux1 = union(Aux1,Aux3);


    % Combined abnormal trials in ascendent order
    % L=0;
    % for j=1:length(Aux1)
    %     L=L+1;
    %     for k=L:length(Aux1)
    %         if Aux1(k)>Aux1(j)
    %             aux=Aux1(j);
    %             Aux1(j)=Aux1(k);
    %             Aux1(k)=aux;
    %         end
    %     end
    % end
    % Aux1 = flip(Aux1);

    Qtarchive = Qt;
    Qparchive = Qp;

    Qt(:,Aux1) = [];
    Qp(:,Aux1) = [];

    figure
    imagesc(Qt')
    set(gca, 'Ydir','normal')
    colorbar

    figure
    imagesc(Qp')
    set(gca, 'Ydir','normal')
    colorbar

    span = 1000;




    % Klist = [4 5 6 7 8 10 12];
    Klist = find(Qp(1,:)>0.5);
    figure(1001)
    hold on
    set(gca,'fontsize',24);
    xlabel('t [sec]');
    ylabel('S');
    figure(1002)
    hold on
    set(gca,'fontsize',24);
    xlabel('t [sec]');
    ylabel('S');
    for k=1:length(Klist)
        [Qtmax,jmax] = max(Qt(:,Klist(k)));
        tc = Qt(1:jmax,Klist(k));
        Pc = Qp(1:jmax,Klist(k));
        figure(1001)
        plot(tc,Pc,'-b','linewidth',2)
        Pcs = smooth(Pc,span);
        figure(1002)
        plot(tc,Pcs,'-r','linewidth',2);
        max(diff(tc))        
    end

    % Function interpolation on a uniform t-axis (dt is different along each
    % t-vector and across trials)
    % dt is made approximately equal to the mean of the dt's

    dt = 0.05;
    ta = 0:dt:600;
    Pa = zeros(1,length(ta));
    Qa = zeros(length(ta),length(Klist));
    figure(1003)
    hold on
    set(gca,'fontsize',24);
    xlabel('t [sec]');
    ylabel('S');
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
        [k,size(Pa)]
        Qa(:,k) = Pa;
        figure(1003)
        plot(ta,Pa,'-g','linewidth',2);      
    end

    figure
    imagesc(Qa')
    set(gca, 'Ydir','normal')
    colorbar

    Pamean = zeros(1,length(ta));
    for j=1:length(ta)
       Pamean(j) = mean(Qa(j,:)); 
    end
    figure(1001)
    plot(ta,Pamean,'r','linewidth',2)

    eta = (Pamean(2)-Pamean(1))/dt;

    % Analytical solution
    jbar = find(isnan(Pamean),1,'first')-1;
    tbar = t(jbar);
    xbar = Pamean(jbar);
    xo = Pamean(1);
    lda1 = 0.045;
    lda2 = 0.005;
    C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
    C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);    
    Pameananl = C1*exp(-lda1*ta)-C2*exp(-lda2*ta)+xbar;
    plot(ta,Pameananl,'-g','linewidth',2)

    % Gradient descents
    xdata = Pamean(1:jbar);
    tdata = ta(1:jbar);

    % Klist = find(Qp(1,:)>1);

    % lda1 = 0.045;
    % lda2 = 0.005;
    % Klist = find(Qp(1,:)>0.5);
    % lda1 = 0.0477
    % lda2 = 0.0063

    % k = 4;
    % Graphsmooth(Qt,Qp,k,span);
    % k = 6;
    % [Qtmax,jmax] = max(Qt(:,Klist(k)));
    % tc = Qt(1:jmax,Klist(k));
    % Pc = Qp(1:jmax,Klist(k));  
    % figure
    % hold on
    % plot(tc,Pc,'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % dt = 0.05;
    % ta = 0:dt:600;
    % Pa = zeros(1);
    % Pa(1) = Pc(1);
    % for j=2:length(ta)
    %     jcant = find(tc<ta(j),1,'last');
    %     tcant = tc(jcant);
    %     jcpost = find(tc>ta(j),1,'first');
    %     tcpost = tc(jcpost);
    %     Pcant = Pc(jcant);
    %     Pcpost = Pc(jcpost);
    %     Pa(j) = Pcant+(ta(j)-tcant)/(tcpost-tcant)*(Pcpost-Pcant);
    % end
    % plot(ta,Pa,'-g','linewidth',2)




    % %%%%%
    % k = 4;
    % [Qtmax,jmax] = max(Qt(:,k));
    % tc = Qt(1:jmax,k);
    % Pc = Qp(1:jmax,k);
    % figure
    % hold on
    % plot(tc,Pc,'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 1000;
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % 
    % k = 5;
    % [Qtmax,jmax] = max(Qt(:,k));
    % tc = Qt(1:jmax,k);
    % Pc = Qp(1:jmax,k);
    % figure
    % hold on
    % plot(tc,Pc,'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 1000;
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);

    % %%%%%
    % k=9;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 1000;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % %%%%%
    % 
    % k=11;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 800;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % %%%%%
    % 
    % k=12;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 1000;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % %%%%%
    % 
    % k = 14;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 800;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % 
    % %%%%%
    % 
    % k = 19;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 800;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % %%%%%
    % 
    % k = 61;
    % figure
    % hold on
    % plot(t(M(k)+1:M(k+1)),P(M(k)+1:M(k+1)),'-b','linewidth',2)
    % set(gca,'fontsize',24);
    % xlabel('t [sec]');
    % ylabel('S');
    % 
    % span = 800;
    % tc = t(M(k)+1:M(k+1));
    % Pc = P(M(k)+1:M(k+1));
    % Pcs = smooth(Pc,span);
    % plot(tc,Pcs,'-r','linewidth',2);
    % 
    % 
    
elseif KSE == 2
    
    load('horacio1.mat')
    Id = ID;
    t = X;
    P = Y;
    
    span = 1000;
    
    % Boundaries between data sets
    
    M = 0;
    for j=1:length(t)-1
        if t(j) > t(j+1)        
            M = [M j];
        end    
    end
    
    % Number of trials (scalar)
    Ntr = length(M)-1;
    
    % Size of each data set  (vector)
    Msze = diff(M);
    
    % Matrix Qt: times
    % Matrix Qp: Data
    
    Qt = zeros(max(Msze),Ntr);
    Qp = zeros(max(Msze),Ntr);
    for k=1:Ntr
        Qt(1:Msze(k),k) = t(M(k)+1:M(k+1));
        Qp(1:Msze(k),k) = P(M(k)+1:M(k+1));
    end

    figure
    imagesc(Qt')
    set(gca, 'Ydir','normal')
    colorbar

    figure
    imagesc(Qp')
    set(gca, 'Ydir','normal')
    colorbar
    
    % Exclusions: short data sets
    Aux1 = find(Msze<6000);
    
    Qtarchive = Qt;
    Qparchive = Qp;

    Qt(:,Aux1) = [];
    Qp(:,Aux1) = [];

    figure
    imagesc(Qt')
    set(gca, 'Ydir','normal')
    colorbar

    figure
    imagesc(Qp')
    set(gca, 'Ydir','normal')
    colorbar
    
    
    
%     figure(1002)
%     hold on
%     set(gca,'fontsize',24);
%     xlabel('t [sec]');
%     ylabel('S');
%     for k=1:length(Klist)        
%         Pcs = smooth(Pc,span);
%         figure(1002)
%         plot(tc,Pcs,'-r','linewidth',2);
%         max(diff(tc))        
%     end
    
    % Function interpolation on a uniform t-axis (dt is different along each
    % t-vector and across trials)
    % dt is made approximately equal to the mean of the dt's
    
    Klist = find(Qp(1,:)>min(Qp(1,:)));
    dt = 0.05;
    ta = 0:dt:600;
    
    Qa = InterpolationUniformDt(Qt,Qp,Klist,ta);

    
    %     figure(1003)
    %     hold on
    %     set(gca,'fontsize',24);
    %     xlabel('t [sec]');
    %     ylabel('S');
    

    % Selecting traces according to the initial value of the curve
    
    Qpomin = 3;
    Qpomax = 6;
    Klist = find(Qp(1,:)>Qpomin & Qp(1,:)<Qpomax);
    
    figure(1001)
    hold on
    set(gca,'fontsize',24);
    xlabel('t [sec]');
    ylabel('S');
    for k=1:length(Klist)        
        figure(1001)
        plot(ta,Qa(:,k),'-b','linewidth',2)        
        max(diff(ta))        
    end
    
    % Average
    
    Pamean = zeros(1,length(ta));
    for j=1:length(ta)
       Pamean(j) = mean(Qa(j,Klist)); 
    end
    figure(1001)
    plot(ta,Pamean,'r','linewidth',2)

    % Smoothing: multiple time scales
    
    k = 20;
    span1 = 1000;
    span2 = 10000;
    tbdry = 150;    
    [Qasmoothcomb] = SmoothingTwoTimeScales(Qa,span1,span2,tbdry,k,dt);

    
    figure
    hold on
    plot(ta,Qa(:,k),'-b','linewidth',2);   
    plot(ta,Qasmoothcomb,'-r','linewidth',2)


    
    
end

% k = 3;
% figure
% hold on
% plot(ta,Qa(:,k),'-b','linewidth',2);
% Qasmooth = smooth(Qa(:,k),1000);
% plot(ta,Qasmooth,'-r','linewidth',2)
% [Qauxmin,jmin] = min(Qasmooth(1:floor(200/dt)));
% Qasmooth1 = smooth(Qasmooth(jmin:end),5000);
% plot(ta(jmin:end),Qasmooth1,'-g','linewidth',2)

% % Function interpolation on a uniform t-axis (dt is different along each
%     % t-vector and across trials)
%     % dt is made approximately equal to the mean of the dt's
% 
%     Klist = find(Qp(1,:)>min(Qp(1,:)));
%     dt = 0.05;
%     ta = 0:dt:600;
%     Pa = zeros(1,length(ta));
%     Qa = zeros(length(ta),length(Klist));
% %     figure(1003)
% %     hold on
% %     set(gca,'fontsize',24);
% %     xlabel('t [sec]');
% %     ylabel('S');
%     for k=1:length(Klist)
%         [Qtmax,jmax] = max(Qt(:,Klist(k)));
%         tc = Qt(1:jmax,Klist(k));
%         Pc = Qp(1:jmax,Klist(k));
%         Pa(1) = Pc(1);
%         for j=2:length(ta)        
%             jcant = find(tc<ta(j),1,'last');
%             tcant = tc(jcant);
%             jcpost = find(tc>ta(j),1,'first');
%             tcpost = tc(jcpost);
%             Pcant = Pc(jcant);
%             Pcpost = Pc(jcpost);
%             Pa(j) = Pcant+(ta(j)-tcant)/(tcpost-tcant)*(Pcpost-Pcant);
%         end
%         [k,size(Pa)]
%         Qa(:,k) = Pa;
% %         figure(1003)
% %         plot(ta,Pa,'-g','linewidth',2);      
%     end
% 

% % Smoothing: multiple time scales
%
% k = 30;
% figure
% hold on
% plot(ta,Qa(:,k),'-b','linewidth',2);
% Qasmooth = smooth(Qa(:,k),1000);
% plot(ta,Qasmooth,'-r','linewidth',2)
% [Qauxmin,jmin] = min(Qasmooth(1:floor(200/dt)));
% Qasmooth1 = smooth(Qasmooth(jmin:end),10000);
% plot(ta(jmin:end),Qasmooth1,'-g','linewidth',2)
% Qasmoothcomb = [Qasmooth(1:jmin-1) Qasmooth1];
% plot(ta,Qasmoothcomb,'-g','linewidth',2)
