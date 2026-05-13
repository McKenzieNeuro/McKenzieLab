
    function draw_occupancy(pts, edges, n_bins, base_col, max_alpha)
    % Count points per voxel
    ix = discretize(pts(:,1), edges{1});
    iy = discretize(pts(:,2), edges{2});
    iz = discretize(pts(:,3), edges{3});
    valid = ~isnan(ix)&~isnan(iy)&~isnan(iz);
    ix=ix(valid); iy=iy(valid); iz=iz(valid);

    counts = accumarray([ix iy iz], 1, [n_bins n_bins n_bins]);
    counts = counts / (max(counts(:))+1e-12);   % normalize 0-1

    % Draw each occupied voxel as a small transparent cube patch
    xc = (edges{1}(1:end-1)+edges{1}(2:end))/2;
    yc = (edges{2}(1:end-1)+edges{2}(2:end))/2;
    zc = (edges{3}(1:end-1)+edges{3}(2:end))/2;
    dx = edges{1}(2)-edges{1}(1);
    dy = edges{2}(2)-edges{2}(1);
    dz = edges{3}(2)-edges{3}(1);

    thresh = 0.05;   % skip very sparse voxels
    for xi=1:n_bins
        for yi=1:n_bins
            for zi=1:n_bins
                occ = counts(xi,yi,zi);
                if occ < thresh, continue; end

                % Cube vertices
                x0=xc(xi)-dx/2; x1=xc(xi)+dx/2;
                y0=yc(yi)-dy/2; y1=yc(yi)+dy/2;
                z0=zc(zi)-dz/2; z1=zc(zi)+dz/2;

                % 6 faces of the cube
                verts = [x0 y0 z0; x1 y0 z0; x1 y1 z0; x0 y1 z0;  % bottom
                    x0 y0 z1; x1 y0 z1; x1 y1 z1; x0 y1 z1]; % top
                faces = [1 2 3 4; 5 6 7 8; 1 2 6 5;
                    2 3 7 6; 3 4 8 7; 4 1 5 8];

                patch('Vertices',verts,'Faces',faces,...
                    'FaceColor',base_col,...
                    'FaceAlpha', min(max_alpha, occ*max_alpha),...
                    'EdgeColor','none')
            end
        end
    end
    end
