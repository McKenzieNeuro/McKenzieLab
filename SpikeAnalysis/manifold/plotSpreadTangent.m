control = [];
stim = [];
for i = 1:5
control = [control  ses(i).tensor.block(1).tangent_rep_spont];
stim = [stim  ses(i).tensor.block(1).tangent_rep_stim];
end
reg = 1e-6;                 % eigenvalue floor for stability (small-N samples)
nS  = numel(control);

logAreaRatio = nan(nS,1);
genEig       = nan(nS,2);   % [lambda1 lambda2], stim/control variance ratios
dropped      = false(nS,1); % samples skipped as degenerate
Zc_all = []; Zs_all = []; Zs_shape = [];

for i = 1:nS
    Xc = control{i};  Xs = stim{i};

    % --- drop degenerate samples: too few points or any NaN coords ---
    if size(Xc,1) < 3 || size(Xs,1) < 3 || ...
       any(isnan(Xc(:)))   || any(isnan(Xs(:)))
        dropped(i) = true;
        continue
    end

    mu_c = mean(Xc,1);
    Sc = cov(Xc);  Ss = cov(Xs);

    W = inv_sqrtm_psd(Sc, reg);     % control whitener, symmetric (ZCA)

    Zc = (Xc - mu_c) * W;           % control  -> ~N(0, I)
    Zs = (Xs - mu_c) * W;           % stim in control-standardized coords
    Zc_all = [Zc_all; Zc];
    Zs_all = [Zs_all; Zs];

    Ss_std = W * Ss * W;            % = Sc^{-1/2} Ss Sc^{-1/2}
    Ss_std = (Ss_std + Ss_std')/2;
    ev = sort(eig(Ss_std),'descend');
    ev = max(real(ev), reg);        % guard tiny-negative / complex from roundoff
    genEig(i,:)      = ev';          % basis-invariant expansion ratios
    logAreaRatio(i)  = 0.5*log(prod(ev));

    % rotation-aligned shape (for pooled shape view; direction is arbitrary)
    [V,~] = eig(Ss_std);  V = fliplr(V);          % major axis first
    Zs_shape = [Zs_shape; (Zs - mean(Zs,1)) * V];
end

% --- report dropped samples ---
if any(dropped)
    fprintf('dropped %d of %d samples (degenerate): %s\n', ...
        nnz(dropped), nS, mat2str(find(dropped)'));
end

keep      = ~dropped;
areaRatio = exp(logAreaRatio(keep));
p = signrank(logAreaRatio(keep));    % paired test vs. no expansion
fprintf('n = %d samples | median area ratio = %.2f  (IQR %.2f-%.2f),  signrank p = %.3g\n', ...
        nnz(keep), median(areaRatio), prctile(areaRatio,25), prctile(areaRatio,75), p);

% ---- Fig 1: pooled standardized scatter ----
figure; hold on; axis equal
scatter(Zc_all(:,1),Zc_all(:,2),4,[.6 .6 .6],'filled','MarkerFaceAlpha',.15);
scatter(Zs_all(:,1),Zs_all(:,2),4,[.85 .2 .2],'filled','MarkerFaceAlpha',.15);
plot_ellipse([0 0], eye(2), 1, 'k-','LineWidth',2);            % control 1-sigma (unit circle)
plot_ellipse(mean(Zs_all,1), cov(Zs_all), 1, 'r-','LineWidth',2);
xlabel('control-standardized axis 1'); ylabel('axis 2');
title('Pooled: control (gray) vs stimulation (red)'); legend off

% ---- Fig 2: per-sample expansion ----
figure;
subplot(1,2,1); histogram(areaRatio,20); xline(1,'k--','LineWidth',1.5);
xlabel('area ratio (stim/control)'); ylabel('# samples'); title('1-\sigma ellipse area ratio');
subplot(1,2,2); 
plot(genEig(:,2),genEig(:,1),'o'); hold on; axis equal
xl=xlim; yl=ylim; lim=[0 max([xl yl 1])]; xlim(lim); ylim(lim);
plot(lim,lim,'k:'); xline(1,'k--'); yline(1,'k--');
xlabel('\lambda_2 (minor)'); ylabel('\lambda_1 (major)'); title('expansion eigenvalues');
% ---- local functions (must be at the very end of the script file) ----


figure; hold on; axis equal
scatter(Zc_all(:,1),Zc_all(:,2),4,[.6 .6 .6],'filled','MarkerFaceAlpha',.15);
scatter(Zs_shape(:,1),Zs_shape(:,2),4,[.85 .2 .2],'filled','MarkerFaceAlpha',.15);
plot_ellipse([0 0], eye(2), 1, 'k-','LineWidth',2);              % control unit circle
plot_ellipse([0 0], cov(Zs_shape), 1, 'r-','LineWidth',2);      % aligned stim
xlabel('aligned axis 1 (major)'); ylabel('aligned axis 2 (minor)');
title('Aligned: control (gray) vs stimulation, major axes co-registered');

% ---- aligned: KDE density contours, control vs stim ----
gx = linspace(-4,4,140); gy = linspace(-3,3,140);
[GX,GY] = meshgrid(gx,gy); pts = [GX(:) GY(:)];

fc = reshape(ksdensity(Zc_all,  pts), size(GX));   % control density
fs = reshape(ksdensity(Zs_shape, pts), size(GX));   % aligned stim density

probs = [0.50 0.90];                                % mass enclosed by each contour
lc = hdr_levels(fc, gx, gy, probs);
ls = hdr_levels(fs, gx, gy, probs);

figure; hold on; axis equal
scatter(Zc_all(:,1),Zc_all(:,2),3,[.6 .6 .6],'filled','MarkerFaceAlpha',.08);
scatter(Zs_shape(:,1),Zs_shape(:,2),3,[.85 .2 .2],'filled','MarkerFaceAlpha',.08);
contour(GX,GY,fc,[lc(1) lc(1)],'LineColor',[.25 .25 .25],'LineWidth',.5);  % control 90%
contour(GX,GY,fc,[lc(2) lc(2)],'LineColor',[.25 .25 .25],'LineWidth',4);  % control 50%
contour(GX,GY,fs,[ls(1) ls(1)],'LineColor',[.85 .2 .2],'LineWidth',.5);    % stim 90%
contour(GX,GY,fs,[ls(2) ls(2)],'LineColor',[.85 .2 .2],'LineWidth',4);    % stim 50%
hk = plot(nan,nan,'-','Color',[.25 .25 .25],'LineWidth',1.5);
hr = plot(nan,nan,'-','Color',[.85 .2 .2],'LineWidth',1.5);
legend([hk hr],{'control','stimulation'},'Location','northeast');
xlabel('aligned axis 1 (major)'); ylabel('aligned axis 2 (minor)');
title('Density contours (50%, 90% mass): control vs stimulation');

% ---- local function (append to end of script file) ----
function lev = hdr_levels(f, gx, gy, probs)
    dA = (gx(2)-gx(1))*(gy(2)-gy(1));
    v  = sort(f(:),'descend');
    cm = cumsum(v)*dA;  cm = cm/cm(end);   % cumulative probability mass
    lev = arrayfun(@(p) v(find(cm>=p,1,'first')), probs);
    lev = sort(lev);                        % contour needs ascending levels
end

function W = inv_sqrtm_psd(S, reg)
    S = (S + S')/2; [V,D] = eig(S);
    d = max(diag(D),0) + reg;
    W = V * diag(1./sqrt(d)) * V';
end

function h = plot_ellipse(mu, S, k, varargin)
    t = linspace(0,2*pi,200); [V,D] = eig((S+S')/2);
    xy = k * V * sqrt(D) * [cos(t); sin(t)];
    h = plot(xy(1,:)+mu(1), xy(2,:)+mu(2), varargin{:});
end