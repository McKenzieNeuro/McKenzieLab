function [estimateLabel] = sm_getSeizurePred_Mattis(fname,rusTree,sz_off,ops)

s = dir(fname);
dur = s.bytes/ops.nCh_featureFile/ops.Fs/2;

ts = sz_off:sz_off+1800;
ts = ceil(ts);
%get prediction
estimateLabel =[];
dat1 =[];
for i = ts



    if i<=dur-ops.durFeat
        tim = i;
        features = ops.features(fname,tim,ops);


        dat1 = [dat1;features];

        if mod(i,100)== 0
            [outpred,conf] = predict(rusTree,dat1);
            estimateLabel = [estimateLabel;outpred conf];
            dat1 =[];
        elseif i > (dur- mod(dur,100))
            [outpred,conf] = predict(rusTree,dat1);
            estimateLabel = [estimateLabel;outpred conf];
            dat1 =[];
        end

    end




end





end