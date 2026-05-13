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
%%
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
      save('R:\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat','ses','-v7.3')
      ssix = ssix+1;
    end
 
end

%%
clear ses
% add time resolved
load('R:\McKenzieLab\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
for i = 1:length(ses)

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
 save('R:\McKenzieLab\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat','ses','-v7.3')
 

%%
six = 1;
clear d1 d_subject m_subject stim_cond
for opsin = 1

    switch opsin
        case 1
            load('R:\McKenzieLab\AAC\AAC_DataForSam\AAC_Final\Arch\manifold.mat')
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

                if ~ismember(c,[3 7])
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
    m_subject(:,:,:,six) = nanmean(m1,4);
    d_subject(:,:,six) = nanmean(d1,3);
    stim_cond(six) = opsin;
    six = six+1;
end
end

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

% After running both analysis loops, correlate at block level
delta_flow_all = [];
JS_stim_all    = [];
gain_all       = [];
rho_all        = [];

for i = 1:length(ses)
    for j = 1:length(ses(i).tensor.block)
        tf = ses(i).tensor.block(j).tensor_flow;
        % mean delta_totalSpread across motifs for this block
        dTS = nanmean(tf.delta.delta_totalSpread);
        if isnan(dTS), continue; end

        delta_flow_all = [delta_flow_all; dTS];
        JS_stim_all    = [JS_stim_all;   JSt_stim_session(i,j)];  % store JSt_stim per block
        gain_all       = [gain_all;       gaintt_session(i,j)];
        rho_all        = [rho_all;        rho_stim_t_session(i,j)];
    end
end

[r_JS,   p_JS]   = corr(delta_flow_all, JS_stim_all,  'Type','Spearman','rows','pairwise');
[r_gain, p_gain] = corr(delta_flow_all, gain_all,     'Type','Spearman','rows','pairwise');
[r_rho,  p_rho]  = corr(delta_flow_all, rho_all,      'Type','Spearman','rows','pairwise');
%%
nCond = size(d_subject,1);
kp = stim_cond==1;
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
opto_resp = [];stim_rip = [];spont_rip =[];
load('R:\AAC\AAC_DataForSam\AAC_Final\ChR2\manifold.mat')
for i = 1:length(ses)
    try
 cd(ses(i).dirN)
    fils = getAllExtFiles(pwd,'mat',0);

    kp = contains(fils,'celltypes') |  contains(fils,'ripples') |  contains(fils,'optoStim')  |  contains(fils,'spikes') | contains(fils,'cell_metrics');
    fils = fils(kp);
    %load all data

    %load data
    clear allcelltypes ripples spikes optoStim cell_metrics
    for j = 1:length(fils)
        load(fils{j})
    end
[binnedPopRipple,bin_times]=populationMatrix(spikes,1,1,500,optoStim.timestamps(:,1),'zscore',true);

cellType = readtable('celltypes.md','FileType','text');
kp = contains(cellType.Var2,'PYR');
opto_resp = [opto_resp;nanmean(binnedPopRipple(kp,:,:),3)];

[binnedPopRipple,bin_times1]=populationMatrix(spikes,.5,.5,100,ses(i).rip_ts,'zscore',true);
stim_rip = [stim_rip;nanmean(binnedPopRipple(kp,:,ses(i).tensor.block(1).stim_ix),3)];
spont_rip = [spont_rip;nanmean(binnedPopRipple(kp,:,ses(i).tensor.block(1).training_ix),3)];

    end
end
%%
figure
[~,ix] = bestmatch([0 .100],bin_times);
u_resp = mean(opto_resp(:,ix(1):ix(2)),2);
[~,b] = sort(u_resp);
imagesc(bin_times,[],opto_resp(b,:))

figure

imagesc(bin_times,[],stim_rip(b,:))

figure

imagesc(bin_times1,[],stim_rip(b,:)-spont_rip(b,:))

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
%%

%% Extended triplet analysis
% Tests whether post-stim decorrelation reflects movement toward stim representation

bins = [0 logspace(log10(1),log10(3600),10)];

clear dc_sp dc_st anchor_to_stim anchor_to_spont ...
      post_stim_to_stim post_spont_to_spont ...
      pre_sim post_sim cond_label

for i = 1:length(ses)

    % Choose representation — swap these to compare PV vs tensor
    simR_pv     = corr(ses(i).binnedPopRipplez');
    simR_tensor = corr(ses(i).tensorModel.C');

    nRip = length(ses(i).rip_ts);
    in   = ismember(1:nRip, ses(i).PV.block(1).stim_ix)';

    % Set diagonal to NaN
    simR_pv(eye(nRip)==1)     = nan;
    simR_tensor(eye(nRip)==1) = nan;

    dt = repmat(ses(i).rip_ts', nRip, 1) - ...
         repmat(ses(i).rip_ts,  1,    nRip);

    % ----------------------------------------------------------
    % TRIPLET EXTRACTION
    % Anchor: ripple j (spontaneous)
    % Next:   ripple j+1 (stim OR spontaneous)
    % After:  ripple j+2 (spontaneous, always)
    % ----------------------------------------------------------
    % Pre-allocate triplet arrays
    cors_pv     = [];
    cors_tensor = [];
    conds       = [];   % 1 = stim-interpolated, 2 = all-spont
    dts         = [];
    anchor_idx  = [];

    for j = 1:nRip-2

        % Anchor must be spontaneous
        if ~in(j), continue; end

        % After must be spontaneous
        if in(j+2), continue; end

        if in(j+1)
            c = 1;   % stim-interpolated triplet
        elseif ~in(j+1)
            c = 2;   % all-spontaneous triplet
        else
            continue
        end

        cors_pv     = [cors_pv;     simR_pv(j,     [j+1 j+2])];
        cors_tensor = [cors_tensor; simR_tensor(j,  [j+1 j+2])];
        conds       = [conds;       c];
        dts         = [dts;         dt(j, [j+1 j+2])];
        anchor_idx  = [anchor_idx;  j];

    end

    % ----------------------------------------------------------
    % ORIGINAL RESULT: decorrelation of j+2 from anchor
    % dc(:,1) = sim(anchor, j+1)
    % dc(:,2) = sim(anchor, j+2)   <- the decorrelation measure
    % ----------------------------------------------------------
    dc_sp(i,:)     = nanmedian(cors_pv(conds==2, :));
    dc_st(i,:)     = nanmedian(cors_pv(conds==1, :));

    % ----------------------------------------------------------
    % NEW: Is the anchor more similar to the intervening stim
    %      ripple than to the subsequent spontaneous ripple?
    %
    % For stim-interpolated triplets:
    %   anchor_to_stim(j)  = sim(anchor, j+1=stim)
    %   post_stim_sim(j)   = sim(anchor, j+2=spont after stim)
    %
    % For control triplets:
    %   anchor_to_spont(j) = sim(anchor, j+1=spont)
    %   post_spont_sim(j)  = sim(anchor, j+2=spont after spont)
    % ----------------------------------------------------------

    % Stim triplets
    kp_st = conds == 1;
    anchor_to_stim(i)      = nanmedian(cors_pv(kp_st, 1));     % sim(anchor, stim)
    post_stim_sim(i)       = nanmedian(cors_pv(kp_st, 2));     % sim(anchor, post-stim spont)

    anchor_to_stim_t(i)    = nanmedian(cors_tensor(kp_st, 1));
    post_stim_sim_t(i)     = nanmedian(cors_tensor(kp_st, 2));

    % Control triplets
    kp_sp = conds == 2;
    anchor_to_spont(i)     = nanmedian(cors_pv(kp_sp, 1));     % sim(anchor, next spont)
    post_spont_sim(i)      = nanmedian(cors_pv(kp_sp, 2));     % sim(anchor, spont+2)

    anchor_to_spont_t(i)   = nanmedian(cors_tensor(kp_sp, 1));
    post_spont_sim_t(i)    = nanmedian(cors_tensor(kp_sp, 2));

    % ----------------------------------------------------------
    % NEW: Direct test of the pull hypothesis
    % 
    % For each stim-interpolated triplet, compute:
    %   sim(j+1=stim, j+2=post-stim spont)
    % and compare to matched control:
    %   sim(j+1=spont, j+2=post-spont spont)
    %
    % If stim pulls j+2 toward the stim representation,
    % sim(stim, post-stim) should EXCEED sim(spont, post-spont)
    % ----------------------------------------------------------

    cors2_pv     = [];
    cors2_tensor = [];
    conds2       = [];
    dts2         = [];

    for j = 1:nRip-2

        if ~in(j), continue; end
        if in(j+2), continue; end

        if in(j+1)
            c = 1;
        elseif ~in(j+1)
            c = 2;
        else
            continue
        end

        % sim between j+1 and j+2 (the key new comparison)
        cors2_pv     = [cors2_pv;     simR_pv(j+1,     j+2)];
        cors2_tensor = [cors2_tensor; simR_tensor(j+1,  j+2)];
        conds2       = [conds2;       c];
        dts2         = [dts2;         dt(j+1, j+2)];

    end

    % sim(j+1, j+2): stim-to-post vs spont-to-post
    next_to_after_stim(i)  = nanmedian(cors2_pv(conds2==1));
    next_to_after_spont(i) = nanmedian(cors2_pv(conds2==2));

    next_to_after_stim_t(i)  = nanmedian(cors2_tensor(conds2==1));
    next_to_after_spont_t(i) = nanmedian(cors2_tensor(conds2==2));

    fprintf('Session %d done\n', i)
end


fprintf('\n========== PV REPRESENTATION ==========\n')
fprintf('Anchor → stim vs anchor → next-spont (control):\n')
[~,p,~,st] = ttest(anchor_to_stim - anchor_to_spont);
d = nanmean(anchor_to_stim - anchor_to_spont) / ...
    nanstd(anchor_to_stim  - anchor_to_spont);
fprintf('  d=%.3f, t=%.3f, p=%.4f\n', d, st.tstat, p)

fprintf('Post-stim spont sim to anchor vs post-spont sim to anchor:\n')
[~,p,~,st] = ttest(post_stim_sim - post_spont_sim);
d = nanmean(post_stim_sim - post_spont_sim) / ...
    nanstd(post_stim_sim  - post_spont_sim);
fprintf('  d=%.3f, t=%.3f, p=%.4f\n', d, st.tstat, p)

fprintf('Stim-to-post vs spont-to-post (pull hypothesis):\n')
[~,p,~,st] = ttest(next_to_after_stim - next_to_after_spont);
d = nanmean(next_to_after_stim - next_to_after_spont) / ...
    nanstd(next_to_after_stim  - next_to_after_spont);
fprintf('  d=%.3f, t=%.3f, p=%.4f\n', d, st.tstat, p)

fprintf('\n========== TENSOR REPRESENTATION ==========\n')
fprintf('Anchor → stim vs anchor → next-spont:\n')
[~,p,~,st] = ttest(anchor_to_stim_t - anchor_to_spont_t);
d = nanmean(anchor_to_stim_t - anchor_to_spont_t) / ...
    nanstd(anchor_to_stim_t  - anchor_to_spont_t);
fprintf('  d=%.3f, t=%.3f, p=%.4f\n', d, st.tstat, p)

fprintf('Stim-to-post vs spont-to-post (pull hypothesis, tensor):\n')
[~,p,~,st] = ttest(next_to_after_stim_t - next_to_after_spont_t);
d = nanmean(next_to_after_stim_t - next_to_after_spont_t) / ...
    nanstd(next_to_after_stim_t   - next_to_after_spont_t);
fprintf('  d=%.3f, t=%.3f, p=%.4f\n', d, st.tstat, p)


figure('Position', [100 100 1300 500]);
cols = [0.5 0.5 0.5; 0.95 0.45 0.1; 0.2 0.5 0.85];

% --- Panel 1: Original decorrelation result ---
subplot(1,4,1)
x = [dc_sp(:,1) dc_st(:,1) dc_sp(:,2) dc_st(:,2)];
positions = [1 1.4 2.2 2.6];
labels    = {'sp→sp','sp→st','sp→sp+1','sp→st+1'};
for k = 1:4
    bar(positions(k), nanmean(x(:,k)), 0.3, ...
        'FaceColor', cols(1+mod(k,2),:), 'EdgeColor','none'); hold on
    errorbar(positions(k), nanmean(x(:,k)), ...
             nanstd(x(:,k))/sqrt(size(x,1)), 'k', 'LineWidth',1.5)
    plot(positions(k)*ones(size(x,1),1), x(:,k), 'o', ...
         'Color',[0.6 0.6 0.6], 'MarkerFaceColor',[0.8 0.8 0.8], ...
         'MarkerSize', 5)
end
set(gca,'XTick', positions, 'XTickLabel', labels, ...
        'XTickLabelRotation', 30, 'FontSize', 10)
ylabel('Median correlation')
title('Decorrelation result')
xline(1.7,'--k'); box off

% --- Panel 2: Anchor similarity to stim vs control next ---
subplot(1,4,2)
dat = [anchor_to_spont(:) anchor_to_stim(:)];
for k = 1:2
    bar(k, nanmean(dat(:,k)), 0.5, ...
        'FaceColor', cols(k+1,:), 'EdgeColor','none'); hold on
    errorbar(k, nanmean(dat(:,k)), ...
             nanstd(dat(:,k))/sqrt(size(dat,1)), 'k','LineWidth',1.5)
    plot(k*ones(size(dat,1),1), dat(:,k), 'o', ...
         'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',5)
end
set(gca,'XTick',1:2,'XTickLabel',{'anchor→spont','anchor→stim'},...
        'XTickLabelRotation',25,'FontSize',10)
ylabel('Median correlation')
title('Anchor similarity to next event')
[~,p] = ttest(dat(:,1)-dat(:,2));
text(1.5, max(nanmean(dat))*1.05, sprintf('p=%.3f',p), ...
     'HorizontalAlignment','center','FontSize',9)
box off

% --- Panel 3: Pull hypothesis — sim(j+1, j+2) ---
subplot(1,4,3)
dat = [next_to_after_spont(:) next_to_after_stim(:)];
for k = 1:2
    bar(k, nanmean(dat(:,k)), 0.5, ...
        'FaceColor', cols(k+1,:), 'EdgeColor','none'); hold on
    errorbar(k, nanmean(dat(:,k)), ...
             nanstd(dat(:,k))/sqrt(size(dat,1)), 'k','LineWidth',1.5)
    plot(k*ones(size(dat,1),1), dat(:,k), 'o', ...
         'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',5)
end
set(gca,'XTick',1:2,'XTickLabel',{'spont→post-spont','stim→post-stim'},...
        'XTickLabelRotation',25,'FontSize',10)
ylabel('Median correlation')
title({'Pull hypothesis:','sim(event_{n+1}, event_{n+2})'})
[~,p] = ttest(dat(:,1)-dat(:,2));
text(1.5, max(nanmean(dat))*1.05, sprintf('p=%.3f',p), ...
     'HorizontalAlignment','center','FontSize',9)
box off

% --- Panel 4: Tensor version of pull hypothesis ---
subplot(1,4,4)
dat = [next_to_after_spont_t(:) next_to_after_stim_t(:)];
for k = 1:2
    bar(k, nanmean(dat(:,k)), 0.5, ...
        'FaceColor', cols(k+1,:), 'EdgeColor','none'); hold on
    errorbar(k, nanmean(dat(:,k)), ...
             nanstd(dat(:,k))/sqrt(size(dat,1)), 'k','LineWidth',1.5)
    plot(k*ones(size(dat,1),1), dat(:,k), 'o', ...
         'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',5)
end
set(gca,'XTick',1:2,'XTickLabel',{'spont→post-spont','stim→post-stim'},...
        'XTickLabelRotation',25,'FontSize',10)
ylabel('Median correlation (tensor)')
title({'Pull hypothesis','(tensor representation)'})
[~,p] = ttest(dat(:,1)-dat(:,2));
text(1.5, max(nanmean(dat))*1.05, sprintf('p=%.3f',p), ...
     'HorizontalAlignment','center','FontSize',9)
box off

sgtitle('Stim-induced decorrelation: pull hypothesis test', ...
        'FontSize', 13, 'FontWeight', 'bold')
%%


%%
close all
figure
semilogx(bins,nanmean(st_v_st),'k','linewidth',3)
%shg
hold on
bins(1) = .5;
semilogx(bins,nanmean(sp_v_sp),'b','linewidth',3)
hold on
semilogx(bins,nanmean(st_v_sp),'r','linewidth',3)

legend({'Spont v Spont','Stim v Spont'})

%%
figure
plotMeanSEM(1:2,dc_st,'r')
plotMeanSEM(1:2,dc_sp,'k')
%%

%% ============================================================
% REPULSION HYPOTHESIS ANALYSIS
% Tests whether stim events push subsequent spontaneous ripples
% away from the anchor along the anchor→stim axis
% ============================================================

clear repulsion_pv repulsion_tensor basin_cross_stim basin_cross_spont

for i = 1:length(ses)

    % Representations
    R_pv     = ses(i).binnedPopRipplez_norm;   % nRip x nNeurons
    R_tensor = ses(i).tensorModel.C;
    Rtz      = zscore(R_tensor,[],1);
    R_tensor = Rtz ./ (vecnorm(Rtz,2,2) + 1e-12);

    nRip   = length(ses(i).rip_ts);
    in     = ismember(1:nRip, ses(i).PV.block(1).stim_ix)';
    dt     = repmat(ses(i).rip_ts',nRip,1) - ...
             repmat(ses(i).rip_ts, 1, nRip);
    clust  = ses(i).PV.block(1).kmean_id(:);

    % ----------------------------------------------------------
    % BUILD TRIPLETS: anchor(spont) → middle → after(spont)
    % ----------------------------------------------------------

    for rep = 1:2

        switch rep
            case 1; R = R_pv;     label = 'pv';
            case 2; R = R_tensor; label = 'tensor';
        end

        proj_along  = [];   % projection of displacement onto anchor→stim axis
        proj_perp   = [];   % perpendicular component magnitude
        cond_vec    = [];   % 1=stim-interpolated, 2=control
        dt_vec      = [];
        basin_cross = [];   % did the post-event ripple cross a cluster boundary?

        for j = 1:nRip-2

            % Anchor must be spontaneous
            if in(j), continue; end
            % Post-event must be spontaneous
            if in(j+2), continue; end
            % Need valid cluster assignments
            if isnan(clust(j)) || isnan(clust(j+2)), continue; end

            if in(j+1)
                c = 1;   % stim interpolated
            elseif ~in(j+1)
                c = 2;   % all spontaneous
            else
                continue
            end

            % ---- vectors ----
            v_anchor = R(j,:)';           % anchor representation
            v_middle = R(j+1,:)';         % middle event (stim or spont)
            v_after  = R(j+2,:)';         % post-event spontaneous ripple

            % ---- anchor → middle axis (unit vector) ----
            axis_am  = v_middle - v_anchor;
            axis_am  = axis_am / (norm(axis_am) + 1e-12);

            % ---- displacement of post-event from anchor ----
            disp_vec = v_after - v_anchor;

            % ---- project displacement onto anchor→middle axis ----
            % Positive = moved toward middle (attraction)
            % Negative = moved away from middle (repulsion)
            proj = dot(disp_vec, axis_am);

            % ---- perpendicular component ----
            perp = norm(disp_vec - proj * axis_am);

            proj_along  = [proj_along;  proj];
            proj_perp   = [proj_perp;   perp];
            cond_vec    = [cond_vec;    c];
            dt_vec      = [dt_vec;      dt(j, j+2)];

            % ---- basin crossing ----
            % Did the post-event ripple leave the anchor cluster?
            basin_cross = [basin_cross; double(clust(j+2) ~= clust(j))];

        end

        % ----------------------------------------------------------
        % REPULSION TEST 1: projection along anchor→middle axis
        % If repulsion: stim triplets have MORE NEGATIVE projection
        % than control triplets
        % ----------------------------------------------------------

        kp_st = cond_vec == 1;
        kp_sp = cond_vec == 2;

        % IRI-matched median
        [mu_st, mu_sp, perm_p] = iri_matched_median(proj_along, ...
                                     kp_st, kp_sp, dt_vec, 500);

        repulsion.(label).proj_stim(i)  = mu_st;
        repulsion.(label).proj_spont(i) = mu_sp;
        repulsion.(label).proj_delta(i) = mu_st - mu_sp;
        repulsion.(label).perm_p(i)     = perm_p;

        % ----------------------------------------------------------
        % REPULSION TEST 2: perpendicular displacement
        % If stim causes repulsion off the anchor→stim axis,
        % perp component should be larger for stim triplets
        % ----------------------------------------------------------

        [mu_st_p, mu_sp_p, perm_p_p] = iri_matched_median(proj_perp, ...
                                           kp_st, kp_sp, dt_vec, 500);

        repulsion.(label).perp_stim(i)  = mu_st_p;
        repulsion.(label).perp_spont(i) = mu_sp_p;
        repulsion.(label).perp_delta(i) = mu_st_p - mu_sp_p;
        repulsion.(label).perp_p(i)     = perm_p_p;

        % ----------------------------------------------------------
        % REPULSION TEST 3: basin crossing probability
        % If stim causes repulsion, post-stim ripples should be
        % more likely to cross into a different cluster
        % ----------------------------------------------------------

        [mu_st_b, mu_sp_b, perm_p_b] = iri_matched_median(basin_cross, ...
                                           kp_st, kp_sp, dt_vec, 500);

        repulsion.(label).basin_stim(i)  = mu_st_b;
        repulsion.(label).basin_spont(i) = mu_sp_b;
        repulsion.(label).basin_delta(i) = mu_st_b - mu_sp_b;
        repulsion.(label).basin_p(i)     = perm_p_b;

        % ----------------------------------------------------------
        % REPULSION TEST 4: distance from stim ripple
        % If stim causes repulsion, post-stim spont should be
        % FURTHER from the stim ripple than post-spont spont
        % is from its middle event — opposite of pull
        % ----------------------------------------------------------

        dist_to_middle = [];
        for jj = 1:length(proj_along)
            % recompute distance post-event to middle
            % (stored implicitly — recompute in separate pass below)
        end

    end

    fprintf('Session %d done\n', i)
end


reps = {'pv','tensor'};

fprintf('\n======= REPULSION HYPOTHESIS =======\n')
for r = 1:2
    label = reps{r};
    fprintf('\n--- %s ---\n', upper(label))

    % Projection along anchor→stim axis
    d = repulsion.(label).proj_delta;
    [~,p,~,st] = ttest(d);
    cohd = nanmean(d)/nanstd(d);
    fprintf('Projection along anchor→stim axis (neg=repulsion):\n')
    fprintf('  mean delta=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(d), cohd, st.tstat, p)

    % Perpendicular displacement
    d = repulsion.(label).perp_delta;
    [~,p,~,st] = ttest(d);
    cohd = nanmean(d)/nanstd(d);
    fprintf('Perpendicular displacement (pos=broader exploration):\n')
    fprintf('  mean delta=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(d), cohd, st.tstat, p)

    % Basin crossing
    d = repulsion.(label).basin_delta;
    [~,p,~,st] = ttest(d);
    cohd = nanmean(d)/nanstd(d);
    fprintf('Basin crossing probability (pos=more transitions):\n')
    fprintf('  mean delta=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(d), cohd, st.tstat, p)

end


figure('Position',[100 100 1200 500]);
cols = [0.85 0.33 0.1; 0.2 0.5 0.85];   % orange=spont, blue=stim

metrics = {'proj','perp','basin'};
ylabels = {'Projection along anchor→stim axis', ...
           'Perpendicular displacement', ...
           'Basin crossing probability'};
titles  = {'Directional repulsion', ...
           'Off-axis displacement', ...
           'Basin crossing'};

for r = 1:2
    label = reps{r};
    for m = 1:3
        subplot(2, 3, (r-1)*3 + m)

        stim_vals  = repulsion.(label).([metrics{m} '_stim'])';
        spont_vals = repulsion.(label).([metrics{m} '_spont'])';

        dat = [spont_vals stim_vals];
        for k = 1:2
            bar(k, nanmean(dat(:,k)), 0.5, ...
                'FaceColor', cols(k,:), 'EdgeColor','none'); hold on
            errorbar(k, nanmean(dat(:,k)), ...
                     nanstd(dat(:,k))/sqrt(sum(~isnan(dat(:,k)))), ...
                     'k','LineWidth',1.5)
            plot(k*ones(size(dat,1),1)+0.05*randn(size(dat,1),1), ...
                 dat(:,k), 'o', 'Color',[0.6 0.6 0.6], ...
                 'MarkerFaceColor',[0.8 0.8 0.8], 'MarkerSize', 5)
        end

        [~,p] = ttest(dat(:,1) - dat(:,2));
        ymax = max(nanmean(dat)) * 1.15;
        text(1.5, ymax, sprintf('p=%.3f',p), ...
             'HorizontalAlignment','center','FontSize',9)

        if m == 1, yline(0,'--k'); end
        set(gca,'XTick',1:2,'XTickLabel',{'Spont','Stim'}, ...
                'FontSize',10)
        ylabel(ylabels{m},'FontSize',9)
        title(sprintf('%s\n(%s)', titles{m}, upper(label)), 'FontSize',10)
        box off
    end
end

sgtitle('Repulsion hypothesis: post-stim manifold displacement', ...
        'FontSize',13,'FontWeight','bold')

%% ----------------------------------------------------------
% HELPER
% ----------------------------------------------------------
