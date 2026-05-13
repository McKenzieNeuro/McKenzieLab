load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
%%
clear stats
for ix = 1:length(ses)
    stim = ismember(1:length(ses(ix).binnedPopRipple),ses(ix).PV.block(1).stim_ix)';
    X = double(ses(ix).binnedPopRipple_tensor);
    A = ses(ix).tensorModel.A;   % 36 × 20

    [N,T,R] = size(X);
    nComp = 6;
    Z = zeros(size(A,2), T, R);  % components × time × trials

    for r = 1:R
        Z(:,:,r) = A' * squeeze(X(:,:,r));
    end
    Z = Z(1:nComp,:,:);

    [stats(ix)] = compute_timeResolved_geometry(Z,ses(ix).PV.block(1).kmean_id,stim);
end

%%
diagnostics_spont = [];
diagnostics_stim = [];


for n = 1:length(stats)
diagnostics_spont_t = [];
diagnostics_stim_t = [];
for cond = 1:9

    switch cond
        case 1
            field = 'alphaMean';
        case 2
            field = 'omegaMean';
        case 3
            field = 'spiralScore';

      
       
        case 4
            field = 'totalSpread';
        case 5
            field = 'rotFracTotal';
        case 6
            field = 'effectiveRank';


            
    end

    diagnostics_spont_t = [diagnostics_spont_t nanmean(stats(n).spont.(field))];
        diagnostics_stim_t = [diagnostics_stim_t nanmean(stats(n).stim.(field))];

end

diagnostics_spont= [diagnostics_spont;diagnostics_spont_t];

diagnostics_stim= [diagnostics_stim;diagnostics_stim_t];
end

%%


diagnostics = [];

for n = 1:length(stats)
diagnostics_t = [];
for cond = 1:3

    switch cond
        case 1
            field = 'delta_rotFracTotal';
        case 2
            field = 'delta_omegaAlphaRatio';
        case 3
            field = 'delta_spiralScore';

      
            
    end

    diagnostics_t = [diagnostics_t nanmedian(stats(n).delta.(field))];
     

end

diagnostics= [diagnostics;diagnostics_t];

end

%%

(diagnostics_spont)


(diagnostics_stim)














