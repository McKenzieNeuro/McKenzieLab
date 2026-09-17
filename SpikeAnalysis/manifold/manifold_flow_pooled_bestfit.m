function manifold_flow_pooled_bestfit(S, varargin)
% MANIFOLD_FLOW_POOLED_BESTFIT  Pool ripple clusters into a common frame by
% finding each cluster's orientation of best fit (groupwise registration on the
% spontaneous flow field), then averaging. Preserves the vector field and the
% anisotropic density that naive pooling washes out.
%
%   pale pastel  : MEAN SPONTANEOUS DENSITY
%   quiver       : MEAN SPONTANEOUS FLOW  (aligned average)
%   heavy contour: STIMULATED DENSITY
%
% Options (name-value)
%   'Sessions'  indices into S                 (default: all non-empty)
%   'TopOnly'   use only D.top_mot clusters      (default false)
%   'MinSp'     min spont ripples per cluster     (default 8)
%   'Scale'     'spread'|'none' per-cluster scalar scale-norm (default 'spread')
%   'GridN'     common-grid resolution            (default 13)
%   'GridLim'   half-extent of common grid (norm. units, default 2.2)
%   'NAngles'   rotation candidates in [0,2pi)     (default 72)
%   'Reflect'   include handedness flips           (default true)
%   'Iters'     groupwise iterations               (default 6)
%   'StimGroups' fcn handle D -> per-ripple label   (default: all stim = 'stim')
%
% Pins the final (otherwise arbitrary) global orientation so the mean stretch
% axis is horizontal.

p = inputParser;
p.addParameter('Sessions',[]); p.addParameter('TopOnly',false);
p.addParameter('MinSp',8);     p.addParameter('Scale','spread');
p.addParameter('GridN',13);    p.addParameter('GridLim',2.2);
p.addParameter('NAngles',72);  p.addParameter('Reflect',true);
p.addParameter('Iters',6);     p.addParameter('StimGroups',[]);
p.addParameter('Mode','overlay');      % 'overlay' | 'difference' | 'massblob'
p.addParameter('DiffType','pdf');      % 'pdf' | 'massrank'
p.addParameter('MassFrac',0.5);        % enclosed-mass fraction for massblob mode
p.addParameter('ShowFlow',false);      % overlay flow quiver in massblob mode
p.addParameter('Diagnostics',true);    % print + plot alignment diagnostics
% precompute params (used only when a raw `ses` struct array is passed)
p.addParameter('nFold',10);    p.addParameter('fold_use',1);
p.addParameter('PreMinSp',5);  p.addParameter('PreMinSt',2);
p.parse(varargin{:}); o = p.Results;

% Accept either the precomputed S cell array OR the raw `ses` struct array.
if isstruct(S)
    fprintf('Raw ses struct array detected; running precompute_sessions...\n');
    S = precompute_sessions(S, o.nFold, o.fold_use, o.PreMinSp, o.PreMinSt);
end

sess = o.Sessions; if isempty(sess), sess = find(~cellfun(@isempty,S(:)'))'; end

% ── 1. gather per-cluster local data in recentered+scaled frame ───────────
C = struct('Fx',{},'Fy',{},'sp',{},'st',{},'stg',{},'B',{});   % per-cluster bundle
for ii = sess(:)'
    D = S{ii}; if isempty(D), continue; end
    clusters = D.mot_valid; if o.TopOnly && isfield(D,'top_mot'), clusters = D.top_mot; end
    if isempty(o.StimGroups), grp = repmat({'stim'},D.R,1);
    else, g = o.StimGroups(D); if ~iscell(g), g = cellstr(string(g)); end; grp = g(:); end

    for cl = clusters(:)'
        kp_sp = D.clust_id==cl & ~D.stim_vec;
        kp_st = D.clust_id==cl &  D.stim_vec;
        usp = D.umap_all(kp_sp,:);  ust = D.umap_all(kp_st,:);
        Rsp = size(usp,1); if Rsp < o.MinSp, continue; end

        c = mean(usp,1);
        if strcmpi(o.Scale,'spread'), s = mean(std(usp-c,0,1))+1e-9; else, s = 1; end
        sp = (usp - c)/s;                    % recentered+scaled spont anchors
        st = (ust - c)/s;                    % recentered+scaled stim points

        % per-ripple UMAP flow vectors at spont anchors (local J, K->UMAP)
        rip = find(kp_sp);
        [Xp,uXp] = build_local(D, rip, usp);
        vel = zeros(Rsp,2);
        for rr = 1:Rsp
            tr  = squeeze(D.Znorm(:,:,rip(rr)));            % K×T
            J   = local_J(Xp, uXp, usp(rr,:));             % 2×K
            vel(rr,:) = (J*(tr(:,end)-tr(:,1)))'/s;        % scaled velocity
        end

        % scattered interpolants for the flow field (NaN outside hull)
        Fx = scatteredInterpolant(sp(:,1),sp(:,2),vel(:,1),'natural','none');
        Fy = scatteredInterpolant(sp(:,1),sp(:,2),vel(:,2),'natural','none');

        % local linear Jacobian of the flow (for rotation/stretch diagnostics)
        Bc = (sp \ vel).';                  % vel ~ sp*Bc'  ->  Bc is 2x2

        k = numel(C)+1;
        C(k).Fx=Fx; C(k).Fy=Fy; C(k).sp=sp; C(k).st=st; C(k).stg=grp(kp_st);
        C(k).B = Bc;
    end
end
nC = numel(C); if nC<2, error('Need >=2 clusters (got %d).',nC); end

% ── 2. common grid + rotation candidates ─────────────────────────────────
gl = linspace(-o.GridLim,o.GridLim,o.GridN);
[GX,GY] = meshgrid(gl,gl); Gp = [GX(:) GY(:)]; nG = size(Gp,1);
ang = linspace(0,2*pi,o.NAngles+1); ang(end)=[];
Rc = cell(0);
for a = ang, Rc{end+1} = [cos(a) -sin(a); sin(a) cos(a)]; end      %#ok
if o.Reflect, F=[1 0;0 -1]; nb=numel(Rc); for q=1:nb, Rc{end+1}=Rc{q}*F; end; end %#ok

    function Fal = sampleAligned(c, M)        % aligned field on Gp: M*f(M' g)
        gp2 = Gp*M;                            % = (M' * Gp')'
        v   = [c.Fx(gp2(:,1),gp2(:,2)), c.Fy(gp2(:,1),gp2(:,2))];
        Fal = v*M.';                           % rotate vectors into frame
    end

% ── 3. groupwise alignment ────────────────────────────────────────────────
T = sampleAligned(C(1), eye(2));               % template init
Mbest = repmat({eye(2)},nC,1);
for it = 1:o.Iters
    for k = 1:nC
        bestE = inf; bestM = eye(2);
        for q = 1:numel(Rc)
            Fal = sampleAligned(C(k), Rc{q});
            ok  = all(~isnan(Fal),2) & all(~isnan(T),2);
            if nnz(ok) < 0.3*nG, continue; end
            e = sum(sum((Fal(ok,:)-T(ok,:)).^2,2));
            if e < bestE, bestE = e; bestM = Rc{q}; end
        end
        Mbest{k} = bestM;
    end
    % recompute template = nanmean of aligned fields
    Acc = zeros(nG,2); Cnt = zeros(nG,1);
    for k = 1:nC
        Fal = sampleAligned(C(k), Mbest{k});
        ok  = all(~isnan(Fal),2);
        Acc(ok,:) = Acc(ok,:) + Fal(ok,:); Cnt(ok) = Cnt(ok)+1;
    end
    T = Acc ./ max(Cnt,1); T(Cnt==0,:) = NaN;
end

% ── 3b. alignment diagnostics (pre-gauge-fix; rotation is gauge-dependent) ─
Mold = Mbest;                                   % alignment rotations (pre Rfix)
okT  = all(~isnan(T),2);  Tnorm = sum(sum(T(okT,:).^2));
resid   = nan(nC,1);  omega_raw = nan(nC,1);  omega_al = nan(nC,1);
stretch = nan(nC,1);  reflected = false(nC,1);
for k = 1:nC
    Fal = sampleAligned(C(k), Mold{k});
    okk = all(~isnan(Fal),2) & okT;
    resid(k) = sum(sum((Fal(okk,:)-T(okk,:)).^2)) / max(Tnorm,1e-12);
    Bc = C(k).B;
    omega_raw(k) = (Bc(2,1)-Bc(1,2))/2;         % rotation in cluster's raw gauge
    Sg2 = (Bc+Bc.')/2; ev = sort(eig(Sg2),'descend');
    stretch(k)   = (ev(1)-ev(2))/2;             % stretch magnitude (gauge-invariant)
    Bal = Mold{k}*Bc*Mold{k}.';                 % flow in the aligned frame
    omega_al(k)  = (Bal(2,1)-Bal(1,2))/2;       % sign flips iff reflected
    reflected(k) = det(Mold{k}) < 0;
end
fracRef = mean(reflected);
sgn = sign(omega_al); sgn(sgn==0)=1; consAl = max(mean(sgn>0),mean(sgn<0));
if o.Diagnostics
    fprintf('\n── alignment diagnostics (%d clusters) ──\n',nC);
    fprintf('  reflected during alignment : %.0f%%\n',100*fracRef);
    fprintf('  |omega| (rotation)  median = %.3g   [gauge-invariant]\n',median(abs(omega_raw)));
    fprintf('  stretch            median = %.3g   [gauge-invariant, alignment uses this]\n',median(stretch));
    if median(abs(omega_raw)) > 2*median(stretch)
        fprintf('  ** WARNING: |omega| >> stretch -> orientation is weakly identified;\n');
        fprintf('     the spiral-vs-node character is largely gauge-ambiguous in UMAP.\n');
    end
    fprintf('  aligned-omega sign consistency = %.0f%% (50%%=no shared handedness)\n',100*consAl);
    fprintf('  post-alignment residual  median = %.2f  (lower = clusters co-orient)\n\n',median(resid));
end

% ── 4. pin global orientation: mean stretch axis -> horizontal ────────────
ok = all(~isnan(T),2);
B  = (Gp(ok,:)\T(ok,:)).';                     % fit linear field T ~ Gp*B'
Sg = (B+B.')/2; [Vv,Dd] = eig(Sg);
[~,mi] = max(diag(Dd)); axv = Vv(:,mi);
phi = atan2(axv(2),axv(1)); Rfix = [cos(-phi) -sin(-phi); sin(-phi) cos(-phi)];
for k=1:nC, Mbest{k} = Rfix*Mbest{k}; end
Tn = T*Rfix.';

% ── 5. average field + pool densities ─────────────────────────────────────
P_sp = []; stimMap = containers.Map('KeyType','char','ValueType','any');
for k = 1:nC
    M = Mbest{k};
    P_sp = [P_sp; C(k).sp*M.'];                                       %#ok
    stf = C(k).st*M.'; g = C(k).stg;
    for u = unique(g(:)')
        key = char(u); sub = stf(strcmp(g,key),:);
        if isempty(sub), continue; end
        if isKey(stimMap,key), stimMap(key)=[stimMap(key);sub]; else, stimMap(key)=sub; end
    end
end

% ── 6. plot ───────────────────────────────────────────────────────────────
L = o.GridLim;
gd = linspace(-L,L,140); [GXd,GYd]=meshgrid(gd,gd);
hsp = bw_silverman(P_sp); dsp_raw = kde_raw(P_sp,GXd,GYd,hsp);
pcol=[0.56 0.49 0.74]; ks=keys(stimMap);
mag = hypot(Tn(:,1),Tn(:,2)); sc = 1.6*mean(diff(gl))/(prctile(mag(ok),90)+1e-9);

if strcmpi(o.Mode,'difference')
    % ---- KDE(stim) - KDE(spont), one panel per stim group -----------------
    sp_pdf = dsp_raw/sum(dsp_raw(:));
    sp_lvl = hdr_level(sp_pdf,0.5);                       % spont 50%-mass level
    spF = sp_pdf; if strcmpi(o.DiffType,'massrank'), spF = massrank(sp_pdf); end
    ng = numel(ks);
    fig = figure('Color','w','Position',[80 80 max(520*ng,560) 580]);
    tl  = tiledlayout(fig,1,ng,'Padding','compact','TileSpacing','compact');
    cmap = diverging_bwr(256); axlast=[];
    for j=1:ng
        Pst = stimMap(ks{j});
        st_raw = kde_raw(Pst,GXd,GYd,bw_silverman(Pst)); st_pdf = st_raw/sum(st_raw(:));
        stF = st_pdf; if strcmpi(o.DiffType,'massrank'), stF = massrank(st_pdf); end
        Dmap = stF - spF;  m = max(abs(Dmap(:)))+1e-12;
        ax = nexttile(tl); hold(ax,'on'); axis(ax,'equal');
        imagesc(ax,gd,gd,Dmap); set(ax,'YDir','normal'); caxis(ax,[-m m]); colormap(ax,cmap);
        contour(ax,gd,gd,sp_pdf,[sp_lvl sp_lvl],'LineColor',[.2 .2 .25],'LineWidth',1.6);
        quiver(ax,Gp(ok,1),Gp(ok,2),Tn(ok,1)*sc,Tn(ok,2)*sc,0,'Color',[.12 .12 .18],...
               'LineWidth',0.9,'MaxHeadSize',2.0,'AutoScale','off');
        plot(ax,0,0,'+','Color',[.1 .1 .1],'MarkerSize',11,'LineWidth',1.6);
        ax.XLim=[-L L]; ax.YLim=[-L L]; ax.XTick=[]; ax.YTick=[]; box(ax,'on');
        title(ax,sprintf('%s  \\minus  spont',ks{j}),'FontSize',11,'FontWeight','bold');
        axlast=ax;
    end
    cb = colorbar(axlast);
    if strcmpi(o.DiffType,'massrank'), cb.Label.String='stim \minus spont (enclosed-mass rank)';
    else, cb.Label.String='stim \minus spont (prob. density)'; end
    title(tl,sprintf('Stim relative to spont   (red = stim-enriched; outline = spont 50%% mass)   %d clusters, %d sessions',...
          nC,numel(sess)),'FontSize',11,'FontWeight','bold');
elseif strcmpi(o.Mode,'massblob')
    % ---- matched enclosed-mass blobs: spont vs stim, with area ratio -------
    frac  = o.MassFrac;
    cellA = (gd(2)-gd(1))^2;
    sp_pdf = dsp_raw/sum(dsp_raw(:));
    sp_lvl = hdr_level(sp_pdf,frac);
    area_sp = nnz(sp_pdf>=sp_lvl)*cellA;  cen_sp = mean(P_sp,1);

    fig = figure('Color','w','Position',[100 100 720 700]);
    ax  = axes(fig,'Color','w','Box','on','XColor',[.5 .5 .5],'YColor',[.5 .5 .5]);
    hold(ax,'on'); axis(ax,'equal'); ax.XLim=[-L L]; ax.YLim=[-L L];
    ax.XTickLabel=''; ax.YTickLabel='';
    xlabel(ax,'aligned dim 1  (mean stretch axis)','Color',[.3 .3 .3]);
    ylabel(ax,'aligned dim 2','Color',[.3 .3 .3]);

    if o.ShowFlow
        mag=hypot(Tn(:,1),Tn(:,2)); sc=1.6*mean(diff(gl))/(prctile(mag(ok),90)+1e-9);
        quiver(ax,Gp(ok,1),Gp(ok,2),Tn(ok,1)*sc,Tn(ok,2)*sc,0,'Color',[.7 .7 .75],...
               'LineWidth',0.8,'MaxHeadSize',2.0,'AutoScale','off');
    end

    % spont filled blob
    hSp=gobjects(1);
    for pol = contour_polys(gd, sp_pdf, sp_lvl)
        hSp=patch(ax,pol{1}(:,1),pol{1}(:,2),pcol,'FaceAlpha',0.35,...
                  'EdgeColor',pcol*0.7,'LineWidth',1.6);
    end

    SCOL=[0.80 .10 .30; 0.05 .55 .60; 0.85 .55 .05; 0.30 .25 .75];
    lh=hSp; lt={sprintf('spont  (%.0f%% mass)',100*frac)};  annot={};
    for j=1:numel(ks)
        Pst=stimMap(ks{j}); if size(Pst,1)<o.MinSp, continue; end
        st_raw=kde_raw(Pst,GXd,GYd,bw_silverman(Pst)); st_pdf=st_raw/sum(st_raw(:));
        st_lvl=hdr_level(st_pdf,frac);
        area_st=nnz(st_pdf>=st_lvl)*cellA; cen_st=mean(Pst,1);
        col=SCOL(mod(j-1,4)+1,:); hb=gobjects(1);
        for pol = contour_polys(gd, st_pdf, st_lvl)
            hb=patch(ax,pol{1}(:,1),pol{1}(:,2),col,'FaceAlpha',0.18,...
                     'EdgeColor',col,'LineWidth',2.4);
        end
        % centroid displacement arrow (spont centroid -> stim centroid)
        if norm(cen_st-cen_sp) > 0.03*L
            quiver(ax,cen_sp(1),cen_sp(2),cen_st(1)-cen_sp(1),cen_st(2)-cen_sp(2),0,...
                'Color',col*0.6,'LineWidth',1.8,'MaxHeadSize',1.2,'AutoScale','off');
        end
        plot(ax,cen_st(1),cen_st(2),'o','MarkerFaceColor',col,'MarkerEdgeColor','w','MarkerSize',7);
        if isgraphics(hb)
            lh(end+1)=hb; %#ok
            lt{end+1}=sprintf('%s  (area %.2f\\times, shift %.2f)',ks{j},area_st/area_sp,norm(cen_st-cen_sp)); %#ok
        end
    end
    plot(ax,cen_sp(1),cen_sp(2),'+','Color',[.2 .2 .2],'MarkerSize',12,'LineWidth',1.8);
    legend(ax,lh,lt,'Location','northwest','Box','on','Color','w','FontSize',9);
    title(ax,sprintf('Matched %.0f%%-mass footprints   (%d clusters, %d sessions)',...
          100*frac,nC,numel(sess)),'FontSize',11,'FontWeight','bold','Color',[.15 .15 .15]);
    hold(ax,'off');
else
    % ---- overlay mode (pale density + flow + stim contours) ---------------
    fig = figure('Color','w','Position',[100 100 720 700]);
    ax  = axes(fig,'Color',[0.985 0.985 0.99],'Box','on','GridAlpha',0.10,...
               'XColor',[0.5 .5 .5],'YColor',[0.5 .5 .5]); hold(ax,'on');
    axis(ax,'equal'); ax.XLim=[-L L]; ax.YLim=[-L L]; ax.XTickLabel=''; ax.YTickLabel='';
    xlabel(ax,'aligned dim 1  (mean stretch axis)','Color',[.3 .3 .3]);
    ylabel(ax,'aligned dim 2','Color',[.3 .3 .3]);

    dsp = dsp_raw/max(dsp_raw(:)); img=ones([size(dsp) 3]);
    for ch=1:3, img(:,:,ch)=1-dsp*(1-pcol(ch))*0.45; end
    ih=image(ax,gd,gd,img); ih.AlphaData=min(dsp*0.6,0.6); uistack(ih,'bottom');

    quiver(ax, Gp(ok,1), Gp(ok,2), Tn(ok,1)*sc, Tn(ok,2)*sc, 0, ...
        'Color',[0.30 0.32 0.40],'LineWidth',1.4,'MaxHeadSize',2.5,'AutoScale','off');

    SCOL=[0.80 .10 .30; 0.05 .55 .60; 0.85 .55 .05; 0.30 .25 .75];
    lh=gobjects(0); lt={};
    for j=1:numel(ks)
        Pst=stimMap(ks{j}); if size(Pst,1)<o.MinSp, continue; end
        dst_raw=kde_raw(Pst,GXd,GYd,bw_silverman(Pst)); dst=dst_raw/max(dst_raw(:));
        col=SCOL(mod(j-1,4)+1,:);
        [~,hc]=contour(ax,gd,gd,dst,[.4 .6 .8],'LineColor',col,'LineWidth',2.2);
        lh(end+1)=hc; lt{end+1}=sprintf('%s density',ks{j}); %#ok
    end
    plot(ax,0,0,'+','Color',[.2 .2 .2],'MarkerSize',12,'LineWidth',1.8);
    hd_sp=patch(ax,'XData',nan,'YData',nan,'FaceColor',pcol*0.4+0.6,'EdgeColor','none');
    hd_fl=plot(ax,nan,nan,'-','Color',[.30 .32 .40],'LineWidth',2);
    legend(ax,[hd_sp hd_fl lh],[{'mean spont density','mean spont flow'} lt],...
        'Location','northwest','Box','on','Color','w','FontSize',9);
    title(ax,sprintf('Best-fit-aligned pooled manifold  (%d clusters, %d sessions)',...
          nC,numel(sess)),'FontSize',11,'FontWeight','bold','Color',[.15 .15 .15]);
    hold(ax,'off');
end

% ── 7. diagnostics figure ─────────────────────────────────────────────────
if o.Diagnostics
    fd = figure('Color','w','Position',[120 120 1080 320]);
    t2 = tiledlayout(fd,1,3,'Padding','compact','TileSpacing','compact');

    a1=nexttile(t2);
    histogram(a1,abs(omega_raw),'FaceColor',[.45 .5 .65],'EdgeColor','w'); hold(a1,'on');
    histogram(a1,stretch,'FaceColor',[.85 .55 .2],'EdgeColor','w','FaceAlpha',0.6);
    xline(a1,median(abs(omega_raw)),'Color',[.45 .5 .65],'LineWidth',1.5);
    xline(a1,median(stretch),'Color',[.85 .55 .2],'LineWidth',1.5);
    legend(a1,{'|\omega| rotation','stretch'},'Box','off','FontSize',8);
    title(a1,'Gauge-invariant magnitudes'); xlabel(a1,'rate'); ylabel(a1,'clusters');

    a2=nexttile(t2);
    histogram(a2,omega_raw,'FaceColor',[.6 .6 .6],'EdgeColor','w'); hold(a2,'on');
    histogram(a2,omega_al,'FaceColor',[.80 .10 .30],'EdgeColor','w','FaceAlpha',0.6);
    xline(a2,0,'k:');
    legend(a2,{'\omega raw gauge','\omega after reflection-fold'},'Box','off','FontSize',8);
    title(a2,sprintf('Rotation sign  (%.0f%% reflected)',100*fracRef));
    xlabel(a2,'signed \omega');

    a3=nexttile(t2);
    histogram(a3,resid,'FaceColor',[.3 .55 .45],'EdgeColor','w');
    xline(a3,median(resid),'Color',[.3 .55 .45],'LineWidth',1.5);
    title(a3,'Per-cluster residual'); xlabel(a3,'aligned field MSE / |T|^2'); ylabel(a3,'clusters');

    title(t2,'Alignment diagnostics  (sign of \omega is gauge-dependent; |\omega| & stretch are not)',...
          'FontSize',10.5,'FontWeight','bold');
end
end

% ══════════════════════════════════════════════════════════════════════════
function [Xp,uXp] = build_local(D, rip, usp)
Xp=[]; uXp=[];
for r=1:numel(rip)
    tr=squeeze(D.Znorm(:,:,rip(r)));
    Xp=[Xp, tr(:,2:end)]; uXp=[uXp; repmat(usp(r,:),D.T-1,1)]; %#ok
end
end
function J = local_J(Xp,uXp,anchor)
K=size(Xp,1); d2=sum((uXp-anchor).^2,2); h=median(sqrt(d2))+1e-9;
w=exp(-d2/(2*h^2)); Xa=[Xp;ones(1,size(Xp,2))]; A=(Xa.*w.')*Xa.';
lam=1e-3*trace(A)/(K+1); Jaug=(uXp.'*(Xa.*w.').')/(A+lam*eye(K+1)); J=Jaug(:,1:K);
end
function h = bw_silverman(P)
n=size(P,1); h=std(P,0,1)*n^(-1/6)*1.06; h(h<1e-6)=0.05;
end
function dens = kde2d(P,GX,GY,h)
dens=zeros(size(GX));
for n=1:size(P,1)
    dens=dens+exp(-0.5*(((GX-P(n,1))/h(1)).^2+((GY-P(n,2))/h(2)).^2));
end
dens=dens/max(dens(:));
end
function dens = kde_raw(P,GX,GY,h)         % un-normalized sum of Gaussians
dens=zeros(size(GX));
for n=1:size(P,1)
    dens=dens+exp(-0.5*(((GX-P(n,1))/h(1)).^2+((GY-P(n,2))/h(2)).^2));
end
end
function lvl = hdr_level(pdf,frac)          % density level enclosing `frac` mass
v=sort(pdf(:),'descend'); cs=cumsum(v);
idx=find(cs>=frac*cs(end),1,'first');
if isempty(idx), lvl=min(v(:)); else, lvl=v(idx); end
end
function r = massrank(pdf)                   % HDR-rank field: mass at density>=here
[sv,ord]=sort(pdf(:),'descend'); cm=cumsum(sv)/sum(sv);
r=zeros(numel(pdf),1); r(ord)=cm; r=reshape(r,size(pdf));
end
function cm = diverging_bwr(n)               % blue (neg) -> white -> red (pos)
if nargin<1, n=256; end
h=floor(n/2);
b=[linspace(0.18,1,h)' linspace(0.32,1,h)' linspace(0.72,1,h)'];
r=[linspace(1,0.80,n-h)' linspace(1,0.12,n-h)' linspace(1,0.22,n-h)'];
cm=[b;r];
end
function P = contour_polys(gd, F, lvl)       % polygons of {F>=lvl} at level lvl
Cm = contourc(gd, gd, F, [lvl lvl]); P = {}; k = 1;
while k < size(Cm,2)
    npts = Cm(2,k);
    P{end+1} = [Cm(1,k+1:k+npts).' Cm(2,k+1:k+npts).']; %#ok
    k = k + npts + 1;
end
end
