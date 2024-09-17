function y = doubleExpFit(t,lda1,lda2,eta,xo,xbar)


C1 = (eta+lda2*(xo-xbar))/(lda2-lda1);
C2 = (eta+lda1*(xo-xbar))/(lda2-lda1);
y = C1*exp(-lda1*t)-C2*exp(-lda2*t)+xbar;

end