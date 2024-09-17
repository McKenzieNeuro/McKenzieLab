function pred = myPred1(x0,y0,x)
sigma0 = 0.2;
kparams0 = [3.5, 6.2];


pred = predict(fitrgp(x0,y0),[x]);

