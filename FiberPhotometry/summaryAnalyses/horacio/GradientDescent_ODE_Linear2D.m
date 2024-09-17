function [lda1,lda2] = GradientDescent_ODE_Linear2D(t,xdata)


    % Gradient descent algorithm to approximate an adaptive pattern
    % (overshoot response to a constant input) by the solution to a linear
    % 2D ODE
    
    % Input: adaptive pattern (x, t)
    
    % Output: eigenvalues
    
    % Solution to a linear ODE system with initial conditions x(0)=y(0)=0
    % eta = (xdata(2)-xdata(1))/dt;
    % C1 = (eta-lda2*xbar)/(lda2-lda1);
    % C2 = (eta-lda1*xbar)/(lda2-lda1);
    % x = C1*exp(-lda1*t)-C2*exp(-lda2*t)+xbar;
    
    % learning rule
    
    lrng = 0.00001;
    
    % Initial values
    
    lda1 = 0.045;
    lda2 = 0.005;
     
    % Data: time information
    
    dt = t(2)-t(1);
    
    % Data: steady-state
    
    xbar = xdata(end);
    
    % Initial estimation
    
    eta = (xdata(2)-xdata(1))/dt;
    xo = xdata(1);
    C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
    C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);   
    xest = C1*exp(-lda1*t)-C2*exp(-lda2*t)+xbar;
    plot(t,xest,'--g','linewidth',2);
    E = sum((xest-xdata).^2)/length(t);
    Eorig = E;
    
    figure(1000)
    hold on
    plot(t,xdata,'-b','linewidth',2);
    plot(t,xest,'--r','linewidth',2);
    axis([0 t(end) -0.1 2]);
    xlabel('t');
    ylabel('x');
    set(gca,'fontsize',24);
    
    % Gradient descent iterations
    
     
    lda1 = 0.045;
    lda2 = 0.005;
    dlda = 0.001;
    for k=1:1000
        lda1p = lda1+dlda;
        lda1m = lda1-dlda;
        lda2p = lda2+dlda;
        lda2m = lda2-dlda;
        xlda1p = (eta+lda2*(xo-xbar))/(lda2-lda1p)*exp(-lda1p*t)-(eta+lda1p*(xo-xbar))/(lda2-lda1p)*exp(-lda2*t)+xbar;
        xlda1m = (eta+lda2*(xo-xbar))/(lda2-lda1m)*exp(-lda1m*t)-(eta+lda1m*(xo-xbar))/(lda2-lda1m)*exp(-lda2*t)+xbar;
        xlda2p = (eta+lda2p*(xo-xbar))/(lda2p-lda1)*exp(-lda1*t)-(eta+lda1*(xo-xbar))/(lda2p-lda1)*exp(-lda2p*t)+xbar;
        xlda2m = (eta+lda2m*(xo-xbar))/(lda2m-lda1)*exp(-lda1*t)-(eta+lda1*(xo-xbar))/(lda2m-lda1)*exp(-lda2m*t)+xbar;
        dxdlda1 = (xlda1p-xlda1m)/(2*dlda);
        dxdlda2 = (xlda2p-xlda2m)/(2*dlda);
        lda1 = lda1-lrng*2*sum((xest-xdata).*dxdlda1)/length(t);
        lda2 = lda2-lrng*2*sum((xest-xdata).*dxdlda2)/length(t);    
        C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
        C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);           
        xest = C1*exp(-lda1*t)-C2*exp(-lda2*t)+xbar;
        E = sum((xest-xdata).^2)/length(t);
        [k lda1 lda2 E]
        if E < 0.001
            break
        end
    end
    
    figure(1000)
    plot(t,xest,'--g','linewidth',2);
