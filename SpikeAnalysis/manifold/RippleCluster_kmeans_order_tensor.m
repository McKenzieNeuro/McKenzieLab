topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\AAC\Acute';
topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\PV';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\CCK';

%topDir = 'R:\AAC\AAC_DataForSam\AAC_Final\Arch';
topDir = 'R:\AAC\AAC_DataForSam\AAC_Final\ChR2';
plotIt = false;
d = dir(topDir);
d = {d(cell2mat({d.isdir})).name};
dirs = d(3:end);

dirs = cellfun(@(a) [topDir filesep a],dirs,'uni',0);
kp = ~contains(dirs,'figures');
dirs = dirs(kp);
%%
%clear ses


optsManifold.k        = 5;
optsManifold.knnType  = 'mutual';    % robust under sparsity
optsManifold.nRefs    = 20;
optsManifold.corrType = 'Spearman';

%ssix = 1;
for i = 13%1:length(dirs)

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

            kNN = 10;

            % ----- reference graph -----
            Dsp = squareform(pdist(X,'cosine'));

            [~,idx] = sort(Dsp,2);
            N = size(X,1);
            nLand = min(500, N);          % number of landmarks
            landmarks = randsample(N, nLand);

            W = sparse(N,N);

            for r = 1:N
                W(r,idx(r,2:kNN+1)) = Dsp(r,idx(r,2:kNN+1));
            end

            W = min(W, W');    % symmetrize

            G = graph(W);
            Dref = distances(G, landmarks);
            Dref(isinf(Dref)) = max(Dref(~isinf(Dref)));


            % helper function
            embed_geodesic = @(P) pdist2(P, X(landmarks,:), 'cosine') + ...
                mean(Dref,2)';

            % spontaneous held-out
            Yg = embed_geodesic(Y);

            % stimulation
            Zg = embed_geodesic(Z);

            % ----- compare geometry -----
            q = linspace(0.1,0.9,9);

            d_ref = quantile(pdist(Dref'), q);
            d_ho  = quantile(pdist(Yg),  q);
            d_st  = quantile(pdist(Zg),  q);

            ses(ssix).(cond).block(b).global_geodesic_ratio_spont = ...
                nanmean(d_ho ./ d_ref);

            ses(ssix).(cond).block(b).global_geodesic_ratio_stim = ...
                nanmean(d_st ./ d_ref);


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
                    [a,b_order ] =max(nanmean(binnedPopRipple5z(:,:,kp_base),3),[],2);
                    in_clust = a > .15;
                    if any(in_clust)
                        %get sequence order
                        [~,b_order_heldout ] =max(nanmean(binnedPopRipple5z(:,:,kp_held_out),3),[],2);
                        [~,b_order_stim ] =max(nanmean(binnedPopRipple5z(:,:,kp_stim),3),[],2);


                        ses(ssix).(cond).block(b).order_corr_stim(ii) = corr(b_order(in_clust),b_order_stim(in_clust),'type','spearman','rows','pairwise');
                        ses(ssix).(cond).block(b).order_corr_spont(ii) = corr(b_order(in_clust),b_order_heldout(in_clust),'type','spearman','rows','pairwise');
                        ses(ssix).(cond).block(b).order_shift_spont(ii) = nanmean(b_order(in_clust) - b_order_heldout(in_clust));
                        ses(ssix).(cond).block(b).order_shift_stim(ii) = nanmean(b_order(in_clust) - b_order_stim(in_clust));

                    end
                    %get held out rate
                    tmp = nanmean(binnedPopRipple5z(:,:,kp_held_out),3);

                    ses(ssix).(cond).block(b).ripple_spont{ii} = tmp;


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


                        % -------------------------------------------------
                        % -------------------------------------------------
                        % GEODESIC DISTORTION (LANDMARK METHOD)
                        % -------------------------------------------------

                        kNN = 10;

                        % ----- reference graph -----
                        Dsp = squareform(pdist(X,'cosine'));

                        [~,idx] = sort(Dsp,2);
                        W = inf(size(Dsp));

                        for r = 1:size(Dsp,1)
                            W(r,idx(r,2:kNN+1)) = Dsp(r,idx(r,2:kNN+1));
                        end

                        G = graph(W,'upper');
                        Dgeo_ref = distances(G);

                        % ----- landmark set -----
                        nLand = min(100, size(X,1));
                        landmarks = randsample(size(X,1), nLand);

                        % reference landmark distances
                        Dref = Dgeo_ref(:,landmarks);

                        % helper function
                        embed_geodesic = @(P) pdist2(P, X(landmarks,:), 'cosine') + ...
                            mean(Dref,1);

                        % spontaneous held-out
                        Yg = embed_geodesic(Y);

                        % stimulation
                        Zg = embed_geodesic(Z);

                        % ----- compare geometry -----
                        q = linspace(0.1,0.9,9);

                        d_ref = quantile(pdist(Dref), q);
                        d_ho  = quantile(pdist(Yg),  q);
                        d_st  = quantile(pdist(Zg),  q);

                        ses(ssix).(cond).block(b).geodesic_ratio_spont(ii) = ...
                            nanmean(d_ho ./ d_ref);

                        ses(ssix).(cond).block(b).geodesic_ratio_stim(ii) = ...
                            nanmean(d_st ./ d_ref);



                    end

                else

                    ses(ssix).(cond).block(b).n_noSupport_spont(ii) = mean(kp_held_out);
                    ses(ssix).(cond).block(b).n_noSupport_stim(ii)= mean(kp_stim);

                end



            end


        end
    end
      save('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat','ses','-v7.3')
      ssix = ssix+1;
    end
 
end

%%
clear ses
% add time resolved
load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
for i = 1:length(ses)
    for j = 1:10
        stim = ismember(1:length(ses(i).binnedPopRipple),ses(i).PV.block(j).stim_ix)';
        X = double(ses(i).binnedPopRipple_tensor);
        A = ses(i).tensorModel.A;   % 36 × 20

        [N,T,R] = size(X);
        nComp = 6;
        Z = zeros(size(A,2), T, R);  % components × time × trials

        for r = 1:R
            Z(:,:,r) = A' * squeeze(X(:,:,r));
        end
        Z = Z(1:nComp,:,:);

        ses(i).tensor.block(j).tensor_flow= compute_timeResolved_geometry(Z,ses(i).PV.block(j).kmean_id,stim);
        

    end
end
 save('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat','ses','-v7.3')
clear ses
% add time resolved
load('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat')
for i = 1:length(ses)
    for j = 1:10
        stim = ismember(1:length(ses(i).binnedPopRipple),ses(i).PV.block(j).stim_ix)';
        X = double(ses(i).binnedPopRipple_tensor);
        A = ses(i).tensorModel.A;   % 36 × 20

        [N,T,R] = size(X);
        nComp = 6;
        Z = zeros(size(A,2), T, R);  % components × time × trials

        for r = 1:R
            Z(:,:,r) = A' * squeeze(X(:,:,r));
        end
        Z = Z(1:nComp,:,:);

        ses(i).tensor.block(j).tensor_flow= compute_timeResolved_geometry(Z,ses(i).PV.block(j).kmean_id,stim);
        

    end
end
 save('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat','ses','-v7.3')


%%
six = 1;
clear d1 d_subject m_subject stim_cond
for opsin = 1:2

    switch opsin
        case 1
            load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
        case 2
            load('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat')
    end

    conds  = [{'Rate'},{'Order'},{'NN'},{'GeodesicRatio'},{'Anisotropic Expansion'},{'Volume'},{'OffManifold'},'FlowSpin'];

%conds  = [{'Rate'},{'Order'},{'NN'},{'GeodesicRatio'},{'Rotation'},{'Anisotropic Expansion'},{'Volume'},{'OffManifold',{'NonLinearDeformation'}}];

nSubjects = length(ses);
nFold = 10;

for s = 1:nSubjects
    for f = 1:nFold

        for c = 1:8
            switch c

                case 1

                    cond =  'rate_';
                case 2

                    cond = 'order_corr_';

                case 3
                    cond = 'nearest_neighbor_global_';



                case 4
                    cond = 'geodesic_ratio_';
                case 10
                    cond = 'tangent_angle_';
                case 5
                    cond = 'jacobian_condition_';
                case 6
                    cond = 'jacobian_volume_';
                case 9
                    cond = 'proc_scale_';


                case 7
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

                if c<8
                stim = ses(s).(embField).block(f).(stimField);
                spont = ses(s).(embField).block(f).(spontField);
                else
                    stim = ses(s).tensor.block(f).tensor_flow.stim.rotFracTotal;
                    spont = ses(s).tensor.block(f).tensor_flow.spont.rotFracTotal;
                    
                end

                mu1 = nanmean(stim);
                mu2 = nanmean(spont);
                m1(c,emb,1,f) = mu1;
                m1(c,emb,2,f) = mu2;
                sd  = nanstd([stim(:); spont(:)]);

                
                if ~ismember(c,[2 4 ])
                    d1(c,emb,f) = (mu1 - mu2) / sd;
                else
                    d1(c,emb,f) = (mu2 - mu1) / sd;

                end
               
            end




        end


                

    end
    m_subject(:,:,:,six) = nanmean(m1,4);
    d_subject(:,:,six) = nanmean(d1,3);
    stim_cond(six) = opsin;
    six = six+1;
end
end

%%
nCond = size(d_subject,1);
kp = stim_cond==2;
clear p
for i = 1:nCond
    for j = 1:2
        [~,p(i,j)] = ttest((squeeze(d_subject(i,j,kp))));
    end
end


figure
[a,b] = sort((nanmean((d_subject(:,2,kp)),3)),'descend');

bar(a)
set(gca,'xtick',1:nCond,'xticklabel',conds(b))
hold on
plot(1:nCond,squeeze(d_subject(b,2,kp)),'o')

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

close all
figure
[a,b] = sort((nanmean((d_subject(:,2,kp)),3)),'descend');

bar(a)
set(gca,'xtick',1:nCond,'xticklabel',conds(b))
hold on
plot(1:nCond,squeeze(d_subject(b,2,kp)),'o')

c = 1:nCond;
hold on

text(c(p(b,2)<.01),1.5*ones(sum(p(b,2)<.01),1),'**')
text(c(p(b,2)<.05),1.5*ones(sum(p(b,2)<.05),1),'*')
ylabel('Cohen''s D')

set(gca,'fontsize',20)
ylim([-.5 1.75])

%%
% plot example


i=1;

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
in = InIntervals(ripples.peaks,optoStim.timestamps);
ses(i).dirN = dirs{i};
ses(i).binnedPopRipple_durPV = .1;
ses(i).binnedPopRipple_nBinTensor = 5;
rip_ts = ripples.peaks;
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

kp_noinf = ~any(isinf(zpop1),2);



    [binnedPopRipple5,bin_times]=populationMatrix(spikes,.05, .05, 5,ripples.timestamps(:,1));
    binnedPopRipple5 = binnedPopRipple5(ispyr,:,:);
 

    Xtensor = round(binnedPopRipple5*.02);
    Xtensor= sptensor(Xtensor);
    clear binnedPopRipple5

    % Fit model
    tensor_opts.K = 20;          % start with 4 templates
    tensor_opts.maxIter = 300;
    tensor_opts.lambdaC = 0.05; % optional sparsity on trial weights
    tensor_opts.verbose = true;
    % model = poissonCP_trialWeights(Xtensor, tensor_opts);
    % C = model.C;
    % 
    % % 1. component-wise standardization
    % Cz = zscore(C,[],1);
    % 
    % % 2. unit-length normalization
    % Cn = Cz ./ vecnorm(Cz,2,2);
    % 
    % 
    % % 3. estimate rate axis in normalized space
    % v_rate = (uRate' * Cn)';
    % v_rate = v_rate / norm(v_rate);
    % 
    % % 4. remove rate direction
    % Cn_ratefree = Cn - (Cn * v_rate) * v_rate';
    % 
    % % 5. renormalize
    % Ztensor = Cn_ratefree ./ vecnorm(Cn_ratefree,2,2);



idxB = find(in& kp_noinf);   % trials to test
idxA = setdiff([1:length(in)]',[idxB;find(~kp_noinf)]);   % manifold-defining trials
rIDXA = idxA(randsample(length(idxA),length(idxA)));
nBlock = 10;
blidx = 1:floor(length(rIDXA)/nBlock):length(rIDXA);
[n,bl] = histc(1:length(rIDXA),blidx);

%%

rawData_cluster = zpop1;

[emb_spont, umap_template] = run_umap( ...
    rawData_cluster(idxA,:), ...
    'save_template_file','umap.mat', ...
    'metric','cosine');


emb_stim = run_umap(rawData_cluster(idxB,:), ...
    'template_file','umap.mat', ...
    'metric','cosine');


embed = nan(size(rawData_cluster,1),2);
embed(idxA,:) = emb_spont;
embed(idxB,:) = emb_stim;

max_k = 30;
avg_silhouette = zeros(max_k-4, 1);  % No silhouette for k=1
for k = 5:max_k
    % Cluster using k-means
    idx = kmeans(emb_spont, k,  'Display', 'off');

    % Compute silhouette values
    s = silhouette(emb_spont, idx);

    % Average silhouette score
    avg_silhouette(k-4) = nanmean(s);
end


[~, best_k] = max(avg_silhouette);
optimal_k = best_k + 4;

ix_train = kmeans(emb_spont,optimal_k, 'Replicates', 5, 'Display', 'off');

Cu = zeros(optimal_k, size(emb_spont,2));
for k = 1:optimal_k
    Cu(k,:) = mean(emb_spont(ix_train==k,:),1);
end

% For training data
ix = nan(length(in),1);
ix(idxA) = ix_train;
ix(idxB) = knnsearch(Cu, emb_stim);

optsManifold.k        = 5;
optsManifold.knnType  = 'mutual';    % robust under sparsity
optsManifold.nRefs    = 20;
optsManifold.corrType = 'Spearman';
manifold = nan(length(in),1);
for b = 1:10
    manA = rIDXA(bl ~=b);
    testA = rIDXA(bl ==b);


    tmp = manifoldAnalysis_ABn(rawData_cluster(manA,:), ...
        rawData_cluster(testA,:), optsManifold);



    manifold(testA,:) = tmp.fracOff;
    tmp = manifoldAnalysis_ABn(rawData_cluster(manA,:), ...
        rawData_cluster(idxB,:), optsManifold);
    manifold(idxB,:) = tmp.fracOff;

end



%%
close all
figure; hold on
col = linspecer(11,'hot');
W =[];
for ii = 1:optimal_k
    kp_held_out = ismember(1:size(zpop1,1),testA)' & ix==ii;
    kp_stim =  ismember(1:size(zpop1,1),idxB)' & ix==ii;
    g(ii) = nanmean(uRatez(ix==ii & in))./nanmean(uRatez(ix==ii & ~in));
    if ~isnan(g(ii))
        
    [~,b] = histc(g(ii),prctile(g,0:10:100));
   
    plot3(embed(ix==ii,1), embed(ix==ii,2),zeros(sum(ix==ii),1),'.','color',col{b})
    hold on
    
    %plot3(embed(kp_held_out,1), embed(kp_held_out,2),manifold(kp_held_out),'.','color',col{b},'markersize',12)
    %plot3(embed(kp_stim,1), embed(kp_stim,2),manifold(kp_stim),'.','color','k','markersize',12)
    X = nanmean(embed(kp_held_out,1));
    Y = nanmean(embed(kp_held_out,2));
    Z = 0;
    W(ii) = nanmedian(manifold(kp_stim)) - nanmean(manifold(kp_held_out));
    quiver3(X,Y,Z,0,0,W(ii),'color',col{b},'linewidth',6)
    end
end
plot3(embed(idxB,1), embed(idxB,2),zeros(length(idxB),1),'.','color','k','markersize',20)

colormap(cell2mat(col))
colorbar
axis off
%%


%%

in = ismember(1:length(ses(1).tensor.block(1).kmean_id),ses(1).tensor.block(1).stim_ix)';
figure
[~,b] = sort(W);
ok = double(ses(1).binnedPopRipple_tensor);
zstack = (ok);

for ii = 1:length(b)
stim = nanmean(zstack(:,:,ix==ii & in),3);
end

%%


close all
figure; hold on


[~,iix] = histc(manifold,0:.05:.8);

col = linspecer(max(ix),'hot');
for ii = 1:max(ix)
    kp_held_out = ~ismember(1:size(zpop1,1),idxB)' & ix==ii;
    kp_stim =  ismember(1:size(zpop1,1),idxB)' & ix==ii;

    plot(embed(kp_held_out,1), embed(kp_held_out,2),'.','color',col{ii})
    hold on
    plot(embed(kp_stim,1), embed(kp_stim,2),'x','color',col{ii})
end

colormap(cell2mat(col))
colorbar

%%
close all
figure
nbin = 150;

xedges = linspace(min(embed(:,1)), max(embed(:,1)), nbin+1);
yedges = linspace(min(embed(:,2)), max(embed(:,2)), nbin+1);
[~,~,~,b] = histcn(embed(idxA,:),xedges,yedges);
ok = accumarray(b,manifold(idxA),[nbin+1 nbin+1],@nanmedian,nan);
k = gaussian2Dfilter([ 100 100],.5);
ok = nanconvn(ok,k,'nanout',true);

[~,~,~,b] = histcn(embed(idxB,:),xedges,yedges);
ok1 = accumarray(b,manifold(idxB),[nbin+1 nbin+1],@nanmedian,nan);
k = gaussian2Dfilter([ 100 100],.25);
ok1 = nanconvn(ok1,k);


diffOK = ok1-ok;

imagesc(xedges,yedges,diffOK)
colormap('bluewhitered')
col = colormap;

col = [.5 .5 .5;col];

diffOK(isnan(diffOK)) = min(diffOK(:))-.00001;
figure
imagesc(xedges,yedges,diffOK')

colormap(col)
hold on
plot(embed(idxB,1),embed(idxB,2),'x','color','k')
%%
%plot rates for each ripple
m =[];rate_stim=[];rate_spont =[];rate_stimz=[];fracOff = [];rate_spont_test = []; N_clust_spont = [];N_clust_stim = [];N_clust_test = [];
for i = 1:length(ses)
    for j = 1:length(ses(i).PV.block)




        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).PV.block(j).manifold_stim,[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
        fracOff = [fracOff;tmp];

        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).uRate(ses(i).PV.block(j).stim_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
        rate_stim = [rate_stim;tmp];
        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).test_ix), ...
            ses(i).uRate(ses(i).PV.block(j).test_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
        rate_spont_test = [rate_spont_test;tmp];



        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).training_ix), ...
            ses(i).uRate(ses(i).PV.block(j).training_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);

        rate_spont = [rate_spont;tmp];

        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).uRate(ses(i).PV.block(j).stim_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);

        N_clust_stim = [N_clust_stim;tmp];

        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).test_ix), ...
            ses(i).uRate(ses(i).PV.block(j).test_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);


        N_clust_test = [N_clust_test;tmp];


        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).training_ix), ...
            ses(i).uRate(ses(i).PV.block(j).training_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);

        N_clust_spont = [N_clust_spont;tmp];



    end
end

gain = (rate_stim -rate_spont_test)./rate_spont_test;
gain(isinf(gain)) = nan;

gain = (rate_stim ./rate_spont_test);
gain(isinf(gain)) = nan;

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
    for j = 1:length(ses(i).PV.block)



         fracOff_stim = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).PV.block(j).manifold_stim,[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
        
         order_corr_stim = ses(i).PV.block(j).order_corr_stim;
          
    
      
        rate_stim = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).uRate(ses(i).PV.block(j).stim_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
      
        rate_spont_test = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).test_ix), ...
            ses(i).uRate(ses(i).PV.block(j).test_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
         
        rate_spont = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).training_ix), ...
            ses(i).uRate(ses(i).PV.block(j).training_ix),[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);

        gaint = (rate_stim-rate_spont_test)./rate_spont;
        gaint(isinf(gaint)) = nan;
         corGain_Ordert(j) = corr(order_corr_stim,rate_spont,'type','spearman','rows','pairwise');
        gaintt(j) = nanmean(gaint);
           gaintstdt(j) = nanstd(gaint);
        N_clust_stim = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).uRate(ses(i).PV.block(j).stim_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);

       
        N_clust_test = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).test_ix), ...
            ses(i).uRate(ses(i).PV.block(j).test_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);

        N_clust_spont = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).training_ix), ...
            ses(i).uRate(ses(i).PV.block(j).training_ix),[ses(i).PV.block(j).optimal_k 1],@numel,nan);
        
        [rho_stim_t(j), pval] = corr(N_clust_stim./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric
        [rho_spont_t(j), pval] = corr(N_clust_test./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric

        p_stim = N_clust_stim / nansum(N_clust_stim);
        p_test = N_clust_test / nansum(N_clust_test);
        p_spont = N_clust_spont / nansum(N_clust_spont);

        valid = (p_stim + p_spont) > 0;

        p_stim = p_stim(valid);
        p_spont1 = p_spont(valid);

         M = 0.5 * (p_stim + p_spont1);

        JSt_stim(j) = 0.5 * nansum(p_stim .* log2(p_stim ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));




        valid = (p_test + p_spont) > 0;
    

        p_test = p_test(valid);
          p_spont1 = p_spont(valid);

          M = 0.5 * (p_test + p_spont1);

        JSt_spont(j) = 0.5 * nansum(p_test .* log2(p_test ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));


       
    end
rho_stim(i) = nanmean(rho_stim_t);
rho_spont(i) = nanmean(rho_spont_t);
gain(i) = nanmean(gaintt);
    gainstd(i)  = nanmean(gaintstdt);
       JS_stim(i) = nanmean(JSt_stim);
    JS_spont(i) = nanmean(JSt_spont);
    corGain_Order(i) = nanmedian(corGain_Ordert);
end

[h,p,~,tstat]  =ttest(JS_stim,JS_spont)
[h,p,~,tstat]  =ttest(rho_spont,rho_stim)

%%


clear gain
%plot rates for each ripple
m =[];rate_stim=[];rate_spont =[];rate_stimz=[];fracOff = [];rate_spont_test = []; N_clust_spont = [];N_clust_stim = [];N_clust_test = [];
for i = 1:length(ses)
    for j = 1:length(ses(i).tensor.block)




      
        rate_stim = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).stim_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).stim_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
      
        rate_spont_test = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).test_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).test_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);
         
        rate_spont = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).training_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).training_ix),[ses(i).tensor.block(j).optimal_k 1],@nanmean,nan);

        gaint = (rate_stim-rate_spont_test)./rate_spont;
        gaint(isinf(gaint)) = nan;

        gaintt(j) = nanmean(gaint);

        gaintstdt(j) = nanstd(gaint);
        N_clust_stim = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).stim_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).stim_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);

       
        N_clust_test = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).test_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).test_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);

        N_clust_spont = accumarray(ses(i).tensor.block(j).kmean_id(ses(i).tensor.block(j).training_ix), ...
            ses(i).uRate(ses(i).tensor.block(j).training_ix),[ses(i).tensor.block(j).optimal_k 1],@numel,nan);

        [rho_stim_t(j), pval] = corr(N_clust_stim./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric
        [rho_spont_t(j), pval] = corr(N_clust_test./N_clust_spont,rate_spont, 'Type', 'Spearman','rows','pairwise'); % non-parametric


        p_stim = N_clust_stim / nansum(N_clust_stim);
        p_test = N_clust_test / nansum(N_clust_test);
        p_spont = N_clust_spont / nansum(N_clust_spont);

        valid = (p_stim + p_spont) > 0;

        p_stim = p_stim(valid);
        p_spont1 = p_spont(valid);

         M = 0.5 * (p_stim + p_spont1);

        JSt_stim(j) = 0.5 * nansum(p_stim .* log2(p_stim ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));




        valid = (p_test + p_spont) > 0;
    

        p_test = p_test(valid);
          p_spont1 = p_spont(valid);

          M = 0.5 * (p_test + p_spont1);

        JSt_spont(j) = 0.5 * nansum(p_test .* log2(p_test ./ M)) + ...
            0.5 * nansum(p_spont1 .* log2(p_spont1 ./ M));

 order_corr_stim = ses(i).tensor.block(j).order_corr_stim;
             corGain_Ordert(j) = corr(order_corr_stim,gaint,'type','spearman','rows','pairwise');
    
       


    end

    JS_stim(i) = nanmean(JSt_stim);
    JS_spont(i) = nanmean(JSt_spont);
rho_stim(i) = nanmean(rho_stim_t);
rho_spont(i) = nanmean(rho_spont_t);
gain(i) = nanmean(gaintt);
gainstd(i)  = nanmean(gaintstdt);
    corGain_Order(i) = nanmean(corGain_Ordert);

end
[h,p,~,tstat]  =ttest(rho_stim(kp),rho_spont(kp))


%%

rate_spont =[];rate_stim=[];rate_spontz =[];rate_stimz=[];fracOff = [];
for i = 1:length(ses)
    for j = 1:length(ses(i).PV.block)
        rate_spont = [rate_spont;ses(i).PV.block(j).rate_spont];
        rate_stim = [rate_stim;ses(i).PV.block(j).rate_stim];
        rate_spontz = [rate_spontz;ses(i).PV.block(j).rate_spontz];
        rate_stimz = [rate_stimz;ses(i).PV.block(j).rate_stimz];

        tmp = accumarray(ses(i).PV.block(j).kmean_id(ses(i).PV.block(j).stim_ix), ...
            ses(i).PV.block(j).manifold_stim,[ses(i).PV.block(j).optimal_k 1],@nanmean,nan);
        fracOff = [fracOff;tmp];
    end
end
   %%


 ripple_order_spont =  [];

 ripple_order_stim =  [];
 for i = 1:length(ses)
tmp = [];tmp1=[];
     for j = 1:10

         tmp = cat(3,tmp,nanmean(ses(i).tensor.block(j).ripple_order_spont,3));
     
              tmp1 = cat(3,tmp1,nanmean(ses(i).tensor.block(j).ripple_order_stim,3));

     end

         ripple_order_spont = cat(3,ripple_order_spont,nanmean(tmp,3));
                  ripple_order_stim = cat(3,ripple_order_stim,nanmean(tmp1,3));

 end

%%

%%
load('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
rate_arch = [];
for i = 1:length(ses)

    for j = 1:10

      rate_arch =   [rate_arch;ses(i).PV.block(j).rate_stim-ses(i).PV.block(j).rate_spont];

    end
end


jacobian_volume_delta_ChR2 =[];
jacobian_volume_delta_arch = [];
for i = 1:length(ses)

    for j = 1:10

      jacobian_volume_delta_arch =   [jacobian_volume_delta_arch;ses(i).PV.block(j).jacobian_volume_stim-ses(i).PV.block(j).jacobian_volume_spont];

    end
end
%%
load('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat')
for i = 1:length(ses)

    for j = 1:10

      jacobian_volume_delta_ChR2 =   [jacobian_volume_delta_ChR2;ses(i).PV.block(j).jacobian_volume_stim-ses(i).PV.block(j).jacobian_volume_spont];

    end
end

rate_ChR2 =[];
for i = 1:length(ses)

    for j = 1:10

      rate_ChR2 =   [rate_ChR2;ses(i).PV.block(j).rate_stim-ses(i).PV.block(j).rate_spont];

    end
end

