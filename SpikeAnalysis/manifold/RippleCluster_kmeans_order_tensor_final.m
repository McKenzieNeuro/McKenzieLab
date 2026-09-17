topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\AAC\Acute';
topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\PV';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\CCK';

topDir = 'R:\AAC\AAC_DataForSam\AAC_Final\Arch';
%topDir = 'R:\AAC\AAC_DataForSam\AAC_Final\ChR2';
plotIt = false;
d = dir(topDir);
d = {d(cell2mat({d.isdir})).name};
dirs = d(3:end);

dirs = cellfun(@(a) [topDir filesep a],dirs,'uni',0);
kp = ~contains(dirs,'figures');
dirs = dirs(kp);

kp = ~contains(dirs,'manifold_analysis');
dirs = dirs(kp);
%%


for rep = 1:10 

clear ses


optsManifold.k        = 5;
optsManifold.knnType  = 'mutual';    % robust under sparsity
optsManifold.nRefs    = 20;
optsManifold.corrType = 'Spearman';

ssix = 1;
for i = 1:length(dirs)

    cd(dirs{i})
    fils = getAllExtFiles(pwd,'mat',0);

    kp = contains(fils,'celltypes') |  contains(fils,'ripples') |  contains(fils,'optoStim')  |  contains(fils,'spikes') | contains(fils,'cell_metrics');
    fils = fils(kp);
    %load all data

    %load data
    clear allcelltypes ripples spikes optoStim cell_metrics
    for j = 1:length(fils)
        load(fils{j})
    end

    if exist('allcelltypes')
        ispyr =cellfun(@(a) isstr(a) && contains(a,'pyr'),allcelltypes);

    else
        ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);
    end

    if sum(ispyr)>20
        rip_ts = ripples.peaks;
        in = InIntervals(rip_ts,optoStim.timestamps);

        %grab just the ripples within an hour of stim

        kp_rip = rip_ts > rip_ts(find(in,1,'first'))-3600 & rip_ts < rip_ts(find(in,1,'last'))+3600;
        rip_ts = rip_ts(kp_rip);
        in = in(kp_rip);

        ses(ssix).dirN = dirs{i};
        ses(ssix).binnedPopRipple_durPV = .1;
        ses(ssix).binnedPopRipple_nBinTensor = 5;

        [binnedPopRipple,bin_times]=populationMatrix(spikes,.05,.05,1,rip_ts);




        binSpk = squeeze(binnedPopRipple)';
        binSpk = binSpk(:,ispyr);
        zpop = zscore(binSpk);
        uRatez = nanmean(zpop(:,:),2);
        uRate = nanmean(binSpk(:,:),2);

        zpop_evt = zpop - mean(zpop,2);
        popNorm = sqrt(sum(zpop_evt.^2,2));
        epsNorm = prctile(popNorm,5);
        zpop1 = zpop_evt ./ max(popNorm, epsNorm);

        opts.use_spont_only = true;
        opts.spont_idx = ~in;
        opts.eps = 1e-6;

        % [zpop1, gain, gain_axis] = sm_popMatNorm(binSpk', opts);

        [binnedPopRipple5z,bin_times]=populationMatrix(spikes,.05,.05,5,rip_ts,'zscore',true);
        binnedPopRipple5z = binnedPopRipple5z(ispyr,:,:);
        [binnedPopRipple5,bin_times]=populationMatrix(spikes,.05, .05,5,rip_ts);
        clear spikes
        binnedPopRipple5 = binnedPopRipple5(ispyr,:,:);
        nBin = ses(ssix).binnedPopRipple_nBinTensor;

        Xtensor = round(binnedPopRipple5*.02);
        Xtensor= sptensor(Xtensor);
        clear binnedPopRipple5

        % Fit model
        tensor_opts.K = 20;          % start with 4 templates
        tensor_opts.maxIter = 300;
        tensor_opts.lambdaC = 0.05; % optional sparsity on trial weights
        tensor_opts.verbose = true;
        model = poissonCP_trialWeights(Xtensor, tensor_opts);
        C = model.C;

        % 1. component-wise standardization
        Cz = zscore(C,[],1);

        % 2. unit-length normalization
        Cn = Cz ./ vecnorm(Cz,2,2);


        % 3. estimate rate axis in normalized space
        v_rate = (uRate' * Cn)';
        v_rate = v_rate / norm(v_rate);

        % 4. remove rate direction
        Cn_ratefree = Cn - (Cn * v_rate) * v_rate';

        % 5. renormalize
        Ztensor = Cn_ratefree ./ vecnorm(Cn_ratefree,2,2);


        %save
        ses(ssix).binnedPopRipplez_norm = zpop1;
        ses(ssix).binnedPopRipplez = zpop;
        ses(ssix).binnedPopRipple = binSpk;
        ses(ssix).uRate = uRate;
        ses(ssix).binnedPopRipple_BinTensorDur = .02;
        ses(ssix).binnedPopRipple_tensor = Xtensor;
        ses(ssix).tensor_opts = tensor_opts;
        ses(ssix).tensorModel = model;
        ses(ssix).Ztensor = Ztensor;

        kp_noinf = ~any(isinf(zpop1),2);


        ses(ssix).rip_ts = rip_ts;
        clear ripples


        % do kmean with cross validation


        max_k = 50;


        idxB = find(in& kp_noinf);   % trials to test
        idxA = setdiff([1:length(in)]',[idxB;find(~kp_noinf)]);   % manifold-defining trials
        rIDXA = idxA(randsample(length(idxA),length(idxA)));
        nBlock = 10;

        blidx = 1:floor(length(rIDXA)/nBlock):length(rIDXA);
        [n,bl] = histc(1:length(rIDXA),blidx);
        ses(ssix).kfold = nBlock;


        for emb = 1:2



            switch emb

                case  1
                    rawData_cluster = zpop1;
                    rawData_manifold = binSpk;
                    cond = 'PV';
                case 2
                    rawData_cluster = (Ztensor);
                    rawData_manifold = C;

                    cond = 'tensor';
            end



            allDist = squareform(pdist(rawData_cluster,'cosine'));



            %define k-fold
            for b = 1:nBlock

                manA = rIDXA(bl ~=b);
                testA = rIDXA(bl ==b);


                ses(ssix).(cond).block(b).training_ix =  manA;
                ses(ssix).(cond).block(b).test_ix =  testA;
                ses(ssix).(cond).block(b).stim_ix =  idxB;

                % loop over directories

                tmp = manifoldAnalysis_ABn(rawData_manifold(manA,:), ...
                    rawData_manifold(idxB,:), optsManifold);
                ses(ssix).(cond).block(b).manifold_stim = tmp.fracOff;

                tmp = manifoldAnalysis_ABn(rawData_manifold(manA,:), ...
                    rawData_manifold(testA,:), optsManifold);
                ses(ssix).(cond).block(b).manifold_spont = tmp.fracOff;



                ses(ssix).(cond).block(b).nearest_neighbor_global_stim = (min(allDist(idxB,manA),[],2));
                ses(ssix).(cond).block(b).nearest_neighbor_global_spont = (min(allDist(testA,manA),[],2));




                varN = num2cell(num2str([1:size(rawData_cluster,2)]'),2);

                [old_embeddinga,embedding]  = run_umap(rawData_cluster(manA,:), 'save_template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false, 'parameter_names',varN,'min_dist',.05,'n_neighbors',10,  'metric','cosine');
                [new_embeddinga] = run_umap(rawData_cluster(testA,:), 'template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false,'parameter_names', varN, 'metric','cosine');
                [new_embeddingb] = run_umap(rawData_cluster(idxB,:), 'template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false,'parameter_names', varN, 'metric','cosine');

                umap_proj = nan(size(rawData_cluster,1),2);
                umap_proj(manA,:) = embedding.embedding;
                ses(ssix).(cond).block(b).umap_proj = umap_proj;

                umap_embedded = nan(size(rawData_cluster,1),2);
                umap_embedded(manA,:) = old_embeddinga;
                umap_embedded(testA,:) = new_embeddinga;
                umap_embedded(idxB,:) = new_embeddingb;


                ses(ssix).(cond).block(b).umap_embedded = umap_embedded;






                avg_silhouette = zeros(max_k-4, 1);  % No silhouette for k=1
                for k = 5:max_k
                    % Cluster using k-means
                    idx = kmeans(umap_embedded(manA,:), k,  'Display', 'off');

                    % Compute silhouette values
                    s = silhouette(umap_embedded(manA,:), idx);

                    % Average silhouette score
                    avg_silhouette(k-4) = nanmean(s);
                end


                [~, best_k] = max(avg_silhouette);
                optimal_k = best_k + 4;
                ses(ssix).(cond).block(b).optimal_k = optimal_k;

                ix_train = kmeans(umap_embedded(manA,:),optimal_k, 'Replicates', 5, 'Display', 'off');

                Cu = zeros(optimal_k, size(umap_embedded,2));
                for k = 1:optimal_k
                    Cu(k,:) = mean(umap_embedded(manA(ix_train==k),:),1);
                end

                % For training data
                ix = nan(size(umap_embedded,1),1);
                ix(manA) = ix_train;
                ix(testA) = knnsearch(Cu, umap_embedded(testA,:));
                ix(idxB)  = knnsearch(Cu, umap_embedded(idxB,:));


                ses(ssix).(cond).block(b).kmean_id = ix;

                %get rates per cluster

                ses(ssix).(cond).block(b).rate_spont  =accumarray(ix(kp_noinf&~in),uRate(kp_noinf&~in),[optimal_k 1],@nanmean,nan);
                ses(ssix).(cond).block(b).rate_spontz =accumarray(ix(kp_noinf&~in),uRatez(kp_noinf&~in),[optimal_k 1],@nanmean,nan);
                ses(ssix).(cond).block(b).rate_stim  =accumarray(ix(kp_noinf&in),uRate(kp_noinf&in),[optimal_k 1],@nanmean,nan);
                ses(ssix).(cond).block(b).rate_stimz =accumarray(ix(kp_noinf&in),uRatez(kp_noinf&in),[optimal_k 1],@nanmean,nan);



                kNN = 10;
                X = rawData_cluster(manA,:);
                Y = rawData_cluster(testA,:);
                Z = rawData_cluster(idxB,:);

                % -------------------------------------------------
                % GEODESIC DISTORTION (LANDMARK METHOD)
                % -------------------------------------------------
               outY = gw_distortion_components(X, Y,'minComp', 30, 'nSub', 50, 'nNull', 200);
               % outZ = gw_distortion_components(X,Z,'minComp', 30, 'nSub', 50, 'nNull', 200);
               % 
               %  ses(ssix).(cond).block(b).global_geodesic_spont = outY;
               % 
               % ses(ssix).(cond).block(b).global_geodesic_stim = outZ;


                ses(ssix).(cond).block(b).p_comp_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).p_comp_stim = nan(optimal_k,1);


                ses(ssix).(cond).block(b).order_corr_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).order_corr_spont = nan(optimal_k,1);


                ses(ssix).(cond).block(b).ripple_order_spont = nan(nBin,nBin,optimal_k);
                ses(ssix).(cond).block(b).ripple_order_stim = nan(nBin,nBin,optimal_k);

                ses(ssix).(cond).block(b).nearest_neighbor_cluster_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).nearest_neighbor_cluster_spont = nan(optimal_k,1);


                ses(ssix).(cond).block(b).n_noSupport_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).n_noSupport_stim = nan(optimal_k,1);


                ses(ssix).(cond).block(b).geodesic_ratio_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).geodesic_ratio_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).tangent_angle_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).tangent_angle_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).jacobian_condition_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).jacobian_condition_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).jacobian_volume_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).jacobian_volume_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_translation_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_translation_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_rotation_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_rotation_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_scale_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_scale_spont = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_residual_stim = nan(optimal_k,1);
                ses(ssix).(cond).block(b).proc_residual_spont = nan(optimal_k,1);



                for ii = 1:optimal_k
                    % define probabliity of new data
                    kp_base = ismember(1:size(rawData_cluster,1),manA)' & ix==ii;
                    kp_held_out = ismember(1:size(rawData_cluster,1),testA)' & ix==ii;
                    kp_stim =  ismember(1:size(rawData_cluster,1),idxB)' & ix==ii;



                    if sum(kp_base)>5
                        mu =  nanmean(  rawData_cluster(kp_base,:),1);
                        sigma =  nancov(  rawData_cluster(kp_base,:));
                        epsilon = 1e-4;  % increase if needed
                        sigma = sigma + epsilon * eye(size(sigma,1));

                        new_data = rawData_cluster(kp_stim,:);
                        old_data = rawData_cluster(kp_held_out,:);
                        tmpp = log(mvnpdf(new_data, mu, sigma));
                        tmpp(isinf(tmpp)) = nan;
                        ses(ssix).(cond).block(b).p_comp_stim(ii)   = nanmean(tmpp);
                        tmpp = log(mvnpdf(old_data, mu, sigma));
                        tmpp(isinf(tmpp)) = nan;
                        ses(ssix).(cond).block(b).p_comp_spont(ii)  = nanmean(tmpp);
                        ses(ssix).(cond).block(b).nearest_neighbor_cluster_stim(ii) = nanmean(min(allDist(kp_stim,kp_base),[],2));
                        ses(ssix).(cond).block(b).nearest_neighbor_cluster_spont(ii) = nanmean(min(allDist(kp_held_out,kp_base),[],2));

                        %get sequence order
                       % [a,b_order ] =max(nanmean(binnedPopRipple5z(:,:,kp_base),3),[],2);
                       % in_clust = a > .15;
                      [b_order, rel, in_clust] = seqOrder(binnedPopRipple5z, kp_base, .15, .3);

                        if any(in_clust)
                            %get sequence order
                            [b_order_heldout] = seqOrder(binnedPopRipple5z, kp_held_out, .15, .3);
                            [b_order_stim] = seqOrder(binnedPopRipple5z, kp_stim, .15, .3);

                            %[~,b_order_heldout ] =max(nanmean(binnedPopRipple5z(:,:,kp_held_out),3),[],2);
                            %[~,b_order_stim ] =max(nanmean(binnedPopRipple5z(:,:,kp_stim),3),[],2);


                            ses(ssix).(cond).block(b).order_corr_stim(ii) = corr(b_order(in_clust),b_order_stim(in_clust),'type','spearman','rows','pairwise');
                            ses(ssix).(cond).block(b).order_corr_spont(ii) = corr(b_order(in_clust),b_order_heldout(in_clust),'type','spearman','rows','pairwise');
                            ses(ssix).(cond).block(b).order_shift_spont(ii) = nanmean(b_order(in_clust) - b_order_heldout(in_clust));
                            ses(ssix).(cond).block(b).order_shift_stim(ii) = nanmean(b_order(in_clust) - b_order_stim(in_clust));

                        end
                        %get held out rate
                        tmp = nanmean(binnedPopRipple5z(:,:,kp_held_out),3);

                        ses(ssix).(cond).block(b).ripple_spont{ii} = tmp;

                        b_order = round(b_order);
                        tmp1 = nan(nBin,nBin);
                        for oo = 1:nBin
                            tmp1(:,oo) = accumarray(b_order(in_clust),tmp(in_clust,oo),[nBin 1],@nanmean,nan);

                        end

                        ses(ssix).(cond).block(b).ripple_order_spont(:,:,ii) = tmp1;


                        tmp = nanmean(binnedPopRipple5z(:,:,kp_stim),3);

                        ses(ssix).(cond).block(b).ripple_stim{ii} = tmp;
                        tmp1 = nan(nBin,nBin);
                        for oo = 1:nBin
                            tmp1(:,oo) = accumarray(b_order(in_clust),tmp(in_clust,oo),[nBin 1],@nanmean,nan);

                        end

                        ses(ssix).(cond).block(b).ripple_order_stim(:,:,ii)= tmp1;
                        % ==============================
                        % MANIFOLD DISTORTION DIAGNOSTICS
                        % ==============================

                        X = rawData_cluster(kp_base,:);
                        Y = rawData_cluster(kp_held_out,:);
                        Z = rawData_cluster(kp_stim,:);

                        if sum(kp_stim)>5 && sum(kp_held_out)>5
                            nPC = min(5, size(X,2));

                            % tangent bases
                            [Ux, Sx] = svd(cov(X), 'econ');
                            [Uy, Sy] = svd(cov(Y), 'econ');
                            [Uz, Sz] = svd(cov(Z), 'econ');


                            Ux = Ux(:,1:nPC);
                            Uy = Uy(:,1:nPC);
                            Uz = Uz(:,1:nPC);

                            % Procrustes on bases (same dimensionality)
                            [d_proc_spont, ~, transform_spont] = procrustes(Ux, Uy, ...
                                'Scaling', true, 'Reflection', false);

                            [d_proc_stim, ~, transform_stim] = procrustes(Ux, Uz, ...
                                'Scaling', true, 'Reflection', false);


                            % -------------------------------------------------
                            % 1. LOCAL PROCRUSTES ALIGNMENT
                            % -------------------------------------------------
                            % Align stim → spontaneous


                            ses(ssix).(cond).block(b).proc_residual_spont(ii) = d_proc_spont;
                            ses(ssix).(cond).block(b).proc_scale_spont(ii) = transform_spont.b;
                            ses(ssix).(cond).block(b).proc_rotation_spont(ii) = ...
                                norm(transform_spont.T - eye(size(transform_spont.T)), 'fro');
                            ses(ssix).(cond).block(b).proc_translation_spont(ii) = ...
                                norm(transform_spont.c(1,:));


                            ses(ssix).(cond).block(b).proc_residual_stim(ii) = d_proc_stim;
                            ses(ssix).(cond).block(b).proc_scale_stim(ii) = transform_stim.b;
                            ses(ssix).(cond).block(b).cluster.proc_rotation_stim(ii) = ...
                                norm(transform_stim.T - eye(size(transform_stim.T)), 'fro');
                            ses(ssix).(cond).block(b).proc_translation_stim(ii) = ...
                                norm(transform_stim.c(1,:));



                            % -------------------------------------------------
                            % center
                            muX = mean(X,1);
                            muY = mean(Y,1);
                            muZ = mean(Z,1);
                            Xc = X - muX;
                            Yc = Y - muY;
                            Zc = Z - muZ;
                            % local tangent bases
                            d = min(5, size(X,2));

                            [Ux,~,~] = svd(cov(Xc),'econ');
                            [Uy,~,~] = svd(cov(Yc),'econ');
                            [Uz,~,~] = svd(cov(Zc),'econ');
                            Ux = Ux(:,1:d);
                            Uy = Uy(:,1:d);
                            Uz = Uz(:,1:d);
                            % project data into tangent coordinates
                            Xp = Xc * Ux;   % [Nx × d]
                            Yp = Yc * Uy;   % [Ny × d]
                            Zp = Zc * Uz;   % [Ny × d]

                            % match sample counts via random pairing
                            n = min(size(Xp,1), size(Yp,1));

                            XpY = Xp(randsample(size(Xp,1), n), :);
                            Yp = Yp(randsample(size(Yp,1), n), :);

                            n = min(size(Xp,1), size(Zp,1));
                            XpZ = Xp(randsample(size(Xp,1), n), :);

                            Zp = Zp(randsample(size(Zp,1), n), :);

                            % ridge regression Jacobian
                            lambda = 1e-3;
                            J = (XpY' * XpY + lambda*eye(d)) \ (XpY' * Yp);

                            % diagnostics
                            s_spont = svd(J);


                            J = (XpZ' * XpZ + lambda*eye(d)) \ (XpZ' * Zp);

                            % diagnostics
                            s_stim = svd(J);


                             ses(ssix).(cond).block(b).tangent_rep_spont{ii} = Xp(:,1:2);
                             ses(ssix).(cond).block(b).tangent_rep_stim{ii}  = Zp(:,1:2);

                            ses(ssix).(cond).block(b).jacobian_singular_values_spont{ii} = s_spont;
                            ses(ssix).(cond).block(b).jacobian_condition_spont(ii) = max(s_spont)/min(s_spont);
                            ses(ssix).(cond).block(b).jacobian_volume_spont(ii) = prod(s_spont);



                            ses(ssix).(cond).block(b).jacobian_singular_values_stim{ii} = s_stim;
                            ses(ssix).(cond).block(b).jacobian_condition_stim(ii) = max(s_stim)/min(s_stim);
                            ses(ssix).(cond).block(b).jacobian_volume_stim(ii) = prod(s_stim);


                            % -------------------------------------------------
                            % 3. TANGENT-SPACE ANGLE CHANGE
                            % -------------------------------------------------
                            % PCA subspace rotation

                            nPC = min(5, size(X,2));

                            [Ux,~,~] = svd(cov(X),'econ');
                            [Uy,~,~] = svd(cov(Y),'econ');
                            [Uz,~,~] = svd(cov(Z),'econ');

                            theta_spont = subspace(Ux(:,1:nPC), Uy(:,1:nPC));
                            theta_stim = subspace(Ux(:,1:nPC), Uz(:,1:nPC));

                            ses(ssix).(cond).block(b).tangent_angle_spont(ii) = theta_spont;
                            ses(ssix).(cond).block(b).tangent_angle_stim(ii) = theta_stim;


                          


                        end

                    else

                        ses(ssix).(cond).block(b).n_noSupport_spont(ii) = mean(kp_held_out);
                        ses(ssix).(cond).block(b).n_noSupport_stim(ii)= mean(kp_stim);

                    end



                end


            end
        end
        ssix = ssix+1;
    end



    % ── Precompute once per session ───────────────────────────────────────
    X  = double(ses(i).binnedPopRipple_tensor);   % N x T x R
    A  = ses(i).tensorModel.A;                    % N x K
    C  = ses(i).tensorModel.C;                    % R x K
    [N,T,R] = size(X);
    nComp = 6;

    % Standardization statistics from trial-weight matrix
    Cz = zscore(C, [], 1);
    Cn = Cz ./ (vecnorm(Cz,2,2)+1e-12);           % R x K

    % Rate axis: per-ripple mean firing rate projected through Cn
    % uRate(r) = mean spike count across neurons for ripple r
    uRate = squeeze(sum(mean(X, 1), 2));    % mean across neurons, sum over time    v_rate = (uRate' * Cn)';                       % K x 1
    v_rate = (uRate' * Cn)';
    v_rate = v_rate / (norm(v_rate) + 1e-12);

    % Raw projection
    Zfull = zeros(size(A,2), T, R);
    for r = 1:R
        Zfull(:,:,r) = A' * squeeze(X(:,:,r));
    end

    % Per-timepoint normalization — same pipeline as Zid
    Cz_mu = mean(Cz,1)';   % K x 1
    Cz_sd = std(Cz,0,1)';  % K x 1
    Znorm = zeros(nComp, T, R);
    for r = 1:R
        for t = 1:T
            z   = Zfull(:,t,r);
            zz  = (z - Cz_mu) ./ (Cz_sd + 1e-12);
            zn  = zz / (norm(zz)+1e-12);
            zrf = zn - (zn' * v_rate) * v_rate;    % orthogonal projection
            zrn = zrf / (norm(zrf)+1e-12);
            Znorm(:,t,r) = zrn(1:nComp);
            %Znorm(:,t,r) = z(1:nComp);
        end
    end

    % ── Per-block analysis ────────────────────────────────────────────────
    for j = 1:10
        stim = ismember(1:R, ses(i).tensor.block(j).stim_ix)';
        ses(i).tensor.block(j).tensor_flow = ...
            compute_timeResolved_geometry(Znorm, ...
            ses(i).tensor.block(j).kmean_id, stim);
    end

end


[out1] = baseline_normalized_geometry(ses);
[out2] = correlate_flow_occupancy(ses);


save(['R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold_analysis\manifold_' num2str(rep) '.mat'],'ses','out1','out2','-v7.3')


end


%%

clear d1 d_subject m_subject stim_cond
for rep = 1:10
six = 1;
   
            load(['R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold_analysis\manifold_' num2str(rep) '.mat'])
        
    conds  = [{'Rate'},{'Order'},{'NN'},{'Anisotropic Expansion'},{'Volume'},{'OffManifold'},'FlowSpin'];

    %conds  = [{'Rate'},{'Order'},{'NN'},{'GeodesicRatio'},{'Rotation'},{'Anisotropic Expansion'},{'Volume'},{'OffManifold',{'NonLinearDeformation'}}];

    nSubjects = length(ses);
    nFold = 10;

    for s = 1:nSubjects
        for f = 1:nFold

            for c = 1:7
                switch c

                    case 1

                        cond =  'rate_';
                    case 2

                        cond = 'order_corr_';

                    case 3
                        cond = 'nearest_neighbor_global_';



                   
                    case 4
                        cond = 'jacobian_condition_';
                    case 5
                        cond = 'jacobian_volume_';
                  

                    case 6
                        cond = 'manifold_';
                end



                stimField = [cond 'stim'];

                spontField = [cond 'spont'];


                for emb = 1:2

                    switch emb
                        case 1
                            embField = 'PV';
                        case 2
                            embField = 'tensor';
                    end

                    if c<7
                        stim = ses(s).(embField).block(f).(stimField);
                        spont = ses(s).(embField).block(f).(spontField);
                        mu1 = nanmean(stim);
                        mu2 = nanmean(spont);
                        m1(c,emb,1,f) = mu1;
                        m1(c,emb,2,f) = mu2;
                        sd  = nanstd([stim(:); spont(:)]);

                    else
                        stim = ses(s).tensor.block(f).tensor_flow.stim.totalSpread;
                        spont = ses(s).tensor.block(f).tensor_flow.spont.totalSpread;

                        mu1 = nanmean(stim);
                        mu2 = nanmean(spont);
                        m1(c,emb,1,f) = mu1;
                        m1(c,emb,2,f) = mu2;
                        sd  = nanstd([stim(:); spont(:)]);

                    end
                    if ismember(c,2)
                        delta = nanmean(spont - stim) ;
                    elseif ~ismember(c,[3 6])
                        delta = nanmean(stim - spont) ;
                    else
                        delta = nanmean(stim) - nanmean(spont) ;
                    end



                    if ~ismember(c,[2 4 ])
                        d1(c,emb,f) = delta/ sd;
                    else
                        d1(c,emb,f) = delta / sd;

                    end

                end




            end




        end
        m_subject(:,:,:,six,rep) = nanmedian(m1,4);
        d_subject(:,:,six,rep) = nanmedian(d1,3);
   
        six = six+1;
    end
end



%%
m_subject = nanmedian(m_subject,5);
d_subject = nanmedian(d_subject,4);
nCond = size(d_subject,1);

clear p
for i = 1:nCond
    for j = 1:2
        [~,p(i,j)] = ttest((squeeze(d_subject(i,j,:))));
    end
end


figure
[a,b] = sort((nanmean((d_subject(:,2,:)),3)),'descend');

bar(a)
set(gca,'xtick',1:nCond,'xticklabel',conds(b))
hold on
plot(1:nCond,squeeze(d_subject(b,2,:)),'o')

c = 1:nCond;
hold on

text(c(p(b,2)<.01),1.5*ones(sum(p(b,2)<.01),1),'**')
text(c(p(b,2)<.05),1.5*ones(sum(p(b,2)<.05),1),'*')
ylabel('Cohen''s D')

set(gca,'fontsize',20)
ylim([-1 1.75])
%%
figure
[a,b] = sort(nanmean(d_subject(:,1,kp),3),'descend');

bar(a)
set(gca,'xtick',1:nCond,'xticklabel',conds(b))
hold on
plot(1:nCond,(squeeze(d_subject(b,1,kp))),'o')

c = 1:nCond;
hold on

text(c(p(b,1)<.01),3.5*ones(sum(p(b,1)<.01),1),'**')
text(c(p(b,1)<.05),3.5*ones(sum(p(b,1)<.05),1),'*')
ylabel('Cohen''s D')
set(gca,'fontsize',20)
ylim([-1 3.75])

%%

clear gain gainstd rho_stim rho_spont corGain_Order
%plot rates for each ripple
m =[];rate_stim=[];rate_spont =[];rate_stimz=[];fracOff = [];rate_spont_test = []; N_clust_spont = [];N_clust_stim = [];N_clust_test = [];
for i = 1:length(ses)
    JSt_stim = [];
    JSt_spont = [];
    gaintt = [];
    gaintstdt =[];
    corGain_Ordert = [];
    rho_stim_t= [];
    rho_spont_t = [];
    for j = 1:length(ses(i).tensor.block)


        %
        % % fracOff_stim = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).stim_ix), ...
        % %     ses(i).tensor.block(j).manifold_stim,[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
        %
        %  %order_corr_stim = ses(i).tensor.block(j).order_corr_stim;
        %
        %
        %
        %  %rate_stim = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).stim_ix), ...
        %      ses(i).uRate(ses(i).tensor.block(j).stim_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
        %
        %  rate_spont_test = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).test_ix), ...
        %      ses(i).uRate(ses(i).tensor.block(j).test_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
        %
        %  rate_spont = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).training_ix), ...
        %      ses(i).uRate(ses(i).tensor.block(j).training_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
        %
        %  gaint = (rate_stim-rate_spont_test)./rate_spont;
        %  gaint(isinf(gaint)) = nan;
        %  corGain_Ordert(j) = corr(order_corr_stim,rate_spont,'type','spearman','rows','pairwise');
        %  gaintt(j) = nanmean(gaint);
        %  gaintstdt(j) = nanstd(gaint);


        N_clust_stim = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).stim_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).stim_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);


        N_clust_test = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).test_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).test_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);

        N_clust_spont = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).training_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).training_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);

        %[rho_stim_t(j), pval] = corr(N_clust_stim./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric
        %[rho_spont_t(j), pval] = corr(N_clust_test./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric

        p_stim = N_clust_stim / nansum1(N_clust_stim);
        p_test = N_clust_test / nansum1(N_clust_test);
        p_spont = N_clust_spont / nansum1(N_clust_spont);

        valid = (p_stim + p_spont) >= 0;

        p_stim = p_stim(valid);
        p_spont1 = p_spont(valid);

        M = 0.5 * (p_stim + p_spont1);

        JSt_stim(j) = 0.5 * nansum(p_stim .* log2(p_stim ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));




        valid = (p_test + p_spont) >= 0;


        p_test = p_test(valid);
        p_spont1 = p_spont(valid);

        M = 0.5 * (p_test + p_spont1);

        JSt_spont(j) = 0.5 * nansum(p_test .* log2(p_test ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));



    end
    %rho_stim(i) = nanmean(rho_stim_t);
    %rho_spont(i) = nanmean(rho_spont_t);
    %gain(i) = nanmean(gaintt);
    %   gainstd(i)  = nanmean(gaintstdt);
    JS_stim(i) = nanmean(JSt_stim);
    JS_spont(i) = nanmean(JSt_spont);
    %  corGain_Order(i) = nanmedian(corGain_Ordert);
end

[h,p,~,tstat]  =ttest(JS_stim,JS_spont)

%[h,p,~,tstat]  =ttest(rho_spont,rho_stim)


%%

mean_spont = [];
mean_stim = [];
delta_rotFracTotal = [];
for i = 1:length(ses)


    for j = 1:10
        delta_rotFracTotal = [delta_rotFracTotal;ses(i).tensor.block(j).tensor_flow.spont.alphaMean];

        mean_spont = cat(3,mean_spont,ses(i).PV.block(j).ripple_order_spont);

        mean_stim =cat(3,mean_stim,ses(i).PV.block(j).ripple_order_stim);

    end
end

%%
bins = [0 logspace(log10(1),log10(3600),10)];
clear st_v_sp sp_v_sp st_v_st dc_sp dc_st spsp
for i = 1:length(ses)
    %simR = (corr(ses(i).tensorModel.C'));
    simR = (corr(ses(i).binnedPopRipplez'));
    nRip = length(ses(i).rip_ts);
    simR(eye(nRip)==1) = nan;

    dt = repmat(ses(i).rip_ts',nRip,1)-repmat(ses(i).rip_ts,1,length(ses(i).rip_ts));
    in = ismember(1:nRip,ses(i).PV.block(1).stim_ix);
    nstim = cumsum(triu(repmat(in,length(in),1),1),2);
    stim_v_spont = simR(in,~in);
    dt_stp = dt(in,~in);

    spont_v_spont = simR(~in,~in);
    dt_sp = dt(~in,~in);
    nstim_sp = nstim(~in,~in);

    stim_v_stim = simR(in,in);
    dt_st = dt(in,in);


    [~,~,~,idx] = histcn([dt_sp(:) nstim_sp(:)],bins,[0 1 inf]);

    kp = all(idx>0,2);
    spont_v_spont = spont_v_spont(:);
    spsp(:,:,i) = accumarray(idx(kp,:),spont_v_spont(kp),[length(bins) 2],@nanmean,nan);
    st_v_st(i,:) = avghist(dt_st(:),stim_v_stim(:),bins,@nanmean);
    st_v_sp(i,:) = avghist(dt_stp(:),stim_v_spont(:),bins,@nanmean);
    sp_v_sp(i,:) = avghist(dt_sp(:),spont_v_spont(:),bins,@nanmean);

    cors = [];conds =[];dts =[];
    for j = 1:nRip-2
        if ~in(j)
            if in(j+1) &  ~in(j+2)
                c=1;
            elseif  ~in(j+1) &  ~in(j+2)

                c=2;

            end
            cors = [cors;simR(j,j+1:j+2)];
            conds = [conds;c];
            dts = [dts;dt(j,j+1:j+2)];
        end

    end

    dc_sp(i,:) = nanmedian(cors(conds==2 ,:));
    dc_st(i,:) = nanmedian(cors(conds==1,:));
    i
end
%

%% ----------------------------------------------------------
% HELPER
% ----------------------------------------------------------



% Fig 1 on manifold

%% fig 2

%% fig 3

%% fig 4
%%
cond = 'order_corr_';
stimField = [cond 'stim'];

                spontField = [cond 'spont'];
     embField = 'tensor';


%get order infor
all_seq_spont =[];
all_seq_stim =[];

for i = 1:length(ses)
    all_seq_spontt = [];
    all_seq_stimt = [];
    for j = 1:10

        all_seq_spontt = cat(3,all_seq_spontt,ses(i).tensor.block(j).ripple_order_spont);
        all_seq_stimt = cat(3,all_seq_stimt,ses(i).tensor.block(j).ripple_order_stim);
    end

    all_seq_spont = cat(3,all_seq_spont,nanmedian(all_seq_spontt,3));
    all_seq_stim = cat(3,all_seq_stim,nanmedian(all_seq_stimt,3));
end


%%

all_corr_stim = [];all_corr_spont = [];
for i = 1:5
for b = 1:10
all_corr_stim = [all_corr_stim;ses(i).(embField).block(b).(stimField)];
all_corr_spont = [all_corr_spont;ses(i).(embField).block(b).(spontField)];


end
end

%%
close all
figure
bar(-1:.1:1,histc(all_corr_spont-all_corr_stim,-1:.1:1))


%%
figure
imagesc(nanmean(all_seq_spont,3),[-.1 .4])
title('spont')
axis off
figure
imagesc(nanmean(all_seq_stim,3),[-.1 .4])
title('stim')
axis off
%%
close all
figure
plotMeanSEM(.02:.02:.1,squeeze(all_seq_spont(1,:,:))','k')

hold on
plotMeanSEM(.02:.02:.1,squeeze(all_seq_stim(1,:,:))','r')
set(gca,'fontsize',16)
ylim([-.2 .6])
figure
plotMeanSEM(.02:.02:.1,squeeze(all_seq_spont(3,:,:))','k')
hold on
plotMeanSEM(.02:.02:.1,squeeze(all_seq_stim(3,:,:))','r')
set(gca,'fontsize',16)
ylim([-.2 .6])
figure
plotMeanSEM(.02:.02:.1,squeeze(all_seq_spont(5,:,:))','k')
hold on
plotMeanSEM(.02:.02:.1,squeeze(all_seq_stim(5,:,:))','r')

set(gca,'fontsize',16)
ylim([-.2 .6])

%%
load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\m218_201106_102900\m218_201106_102900.spikes.cellinfo.mat')
load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\m218_201106_102900\m218_201106_102900.cell_metrics.cellinfo.mat')

if exist('allcelltypes')
    ispyr =cellfun(@(a) isstr(a) && contains(a,'pyr'),allcelltypes);

else
    ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);
end
[binnedPopRipple5z,bin_times]=populationMatrix(spikes,.05,.05,5,ses(3).rip_ts,'zscore',true);

binnedPopRipple5z = binnedPopRipple5z(ispyr,:,:);

%%
for bl = 1
    for c = 7%1:max( ses(3).tensor.block(bl).kmean_id)
close all
        ids = ses(3).tensor.block(bl).kmean_id==c;

        rip_ts  = ses(3).rip_ts(ids);

        stims = ses(3).tensor.block(bl).stim_ix;





        %get sequence order
        [a,b_order ] =max(nanmean(binnedPopRipple5z(:,:,ids),3),[],2);
        in_clust = a > .15;
        if any(in_clust)
            %get sequence order
           
            [~,b_order_stim ] =max(nanmean(binnedPopRipple5z(:,:,stims),3),[],2);

        end

        [~,b ] = sort(b_order);


        figure
        stim_idx = (ismember(1:length(ids),stims));
        stim4 =  ses(3).rip_ts(stim_idx(:) &ids );
        
        for i =[2 3 5 22]

            spks =    cellfun(@(a) a(a>stim4(i)-.1 & a<stim4(i)+.1) -stim4(i),spikes.times,'uni',0 );
            spks = spks(b);
            PipeRaster(spks(in_clust))
            title(['block: ' num2str(bl) ' cluster: ' num2str(c)])
saveas(gcf,['E:\Dropbox\UNM\Papers\AAC\Resubmission\Figures\stim_raster_block' num2str(bl) '_cluster' num2str(c) '_trial' num2str(i) '.eps'],'epsc')
close all
        end




        figure
        for i = [  44 121 186 215]

            spks =    cellfun(@(a) a(a>rip_ts(i)-.1 & a<rip_ts(i)+.1) -rip_ts(i),spikes.times,'uni',0 );
            spks = spks(b);
            PipeRaster(spks(in_clust))
saveas(gcf,['E:\Dropbox\UNM\Papers\AAC\Resubmission\Figures\spont_raster_block' num2str(bl) '_cluster' num2str(c) '_trial' num2str(i) '.eps'],'epsc')
close all
        end

    end
end

%%
figure
plot(flipud(squeeze(m_subject(2,1,:,:))))
figure
hist(squeeze(d_subject(2,2,:)),10)
