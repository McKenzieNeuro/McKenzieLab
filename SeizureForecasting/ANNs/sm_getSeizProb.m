function Prob = sm_getSeizProb(data,model,ops)



feat = ops.feature_fun(data,ops);

Prob = predict(model,feat);




end