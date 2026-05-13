
    function draw_density_cloud2(pts, base_col, max_alpha, marker_scale)
        if size(pts,1) < 3, return; end

        % Grid
        n_grid = 8;
        g = cell(1,3);
        for d=1:size(pts,2)
            lo=prctile(pts(:,d),2); hi=prctile(pts(:,d),98);
            pad=(hi-lo)*0.2; if pad<1e-6, pad=0.05; end
            g{d}=linspace(lo-pad, hi+pad, n_grid);
        end
        [GX,GY] = meshgrid(g{1},g{2});
        grid_pts = [GX(:) GY(:) ];

        % Bandwidth: Silverman's rule per dimension
        n  = size(pts,1);
        bw = std(pts) * n^(-1/5) * 1.5;
        bw(bw<1e-6) = 0.05;

        % Evaluate KDE at each grid point
        dens = zeros(size(grid_pts,1),1);
        for j = 1:n
            d2 = sum(((grid_pts - pts(j,:)) ./ bw).^2, 2);
            dens = dens + exp(-0.5*d2);
        end
        dens = dens / max(dens);

        % Threshold — only show above 15% of peak
        thresh = 0.001;
        kp     = dens > thresh;
        if sum(kp) < 1, return; end

        gp  = grid_pts(kp,:);
        den = dens(kp);

        % Marker size and alpha scale with density
        sz  = 10 * den.^0.5;
        col = repmat(base_col, sum(kp), 1);

        % Draw in alpha-sorted order (low density first, high on top)
        [den_sort, si] = sort(den);
        gp  = gp(si,:);
        sz  = sz(si);

        % scatter3 doesn't support per-point alpha in MATLAB
        % Workaround: draw in batches by density level
        n_levels = 6;
        edges_d  = linspace(thresh, 1.0, n_levels+1);
        for lv = 1:n_levels
            in_lv = den_sort >= edges_d(lv) & den_sort < edges_d(lv+1);
            if sum(in_lv)==0, continue; end
            alpha_lv = max_alpha * ((edges_d(lv)+edges_d(lv+1))/2);
            scatter(gp(in_lv,1), gp(in_lv,2), ...
    sz(in_lv)*500, col(in_lv,:), 'filled', ...
    'MarkerFaceAlpha', alpha_lv*marker_scale, ...
    'MarkerEdgeColor', 'none')
        end
    end
