function R = gto2_run_all(varargin)
% GTO2_RUN_ALL  Regenerate every supplement table, in cascade order, with logging.
%
%   R = gto2_run_all()                        full run
%   R = gto2_run_all('quick', true)           tiny parameters: checks every call runs end to end
%   R = gto2_run_all('only', {'S8','S16'})    run a subset by label (see list below)
%   R = gto2_run_all('outdir', 'myresults')   choose the output folder
%
% Each step is wrapped in try/catch, so one failure does not stop the run. Output
% goes to <outdir>/gto2_run_<timestamp>.log (diary of everything printed) and
% <outdir>/gto2_run_<timestamp>.mat (struct R, saved after every step). R.status
% lists each step with its elapsed time and any error message.
%
% LABELS
%   C       gto2_checks('all')                  invariants C1-C5 (C5 prints Table S1)
%   S2      gto2_gates01('gate0')               sheet proximity
%   S3      gto2_gates01('gate1')               dimension excess, sigma sweep, 3 seeds
%   S4      gto2_paramsweep('rotation')         Gate 2 rotation, 5 seeds
%   S5      gto2_covroute('gradient')           covariance route, curvature gradient
%   S6      gto2_covroute('boundary')           covariance route, boundary, k sweep
%   S7      gto2_curvlaw                        Props 2-3 on quadric patches
%   S8      gto2_manifolds('torus') x nRef      torus, convergence in nRef
%   S9      gto2_paramsweep('anchor')           anchor bias vs free parameters
%   S10     gto2_aniso('anchor')                anisotropy vs the anchor (k = 40)
%   S11     gto2_aniso('fractions')             anisotropy vs G,T,O (k = 100)
%   S12     gto2_aniso('frame')                 anisotropy vs frame and tau^2 (k = 40)
%   S13     gto2_anchorbias                     anchor offset, pure on-manifold probes
%   S14     gto2_manifolds('clifford')          Gate 5, |c| < 1
%   S15     gto2_manifolds('ellipse')           Gate 5 with varying II
%   S16     gto2_paramsweep('cor3')             Eq. 24 vs exact form vs anchor model
%   S17     gto2_offmanifold                    off-manifold recovery, full grid
%   S17fine gto2_offmanifold (suppressive, k = 160, nu = 0.30, O* = 0:0.05:1)

ip = inputParser;
ip.addParameter('quick', false);
ip.addParameter('only', {});
ip.addParameter('outdir', 'gto2_results');
ip.parse(varargin{:});
o = ip.Results;
if ischar(o.only), o.only = {o.only}; end

if exist('OCTAVE_VERSION', 'builtin'), pkg load statistics; end
if ~exist(o.outdir, 'dir'), mkdir(o.outdir); end
stamp   = datestr(now, 'yyyymmdd_HHMMSS');
tag     = ''; if o.quick, tag = '_quick'; end
logfile = fullfile(o.outdir, ['gto2_run_' stamp tag '.log']);
matfile = fullfile(o.outdir, ['gto2_run_' stamp tag '.mat']);
diary(logfile);  diary on;
cleanupObj = onCleanup(@() diary('off'));

Q = o.quick;
steps = {
  'C',       'invariant checks C1-C5 (C5 = Table S1)', ...
      @() gto2_checks('all');
  'S2',      'Gate 0: sheet proximity', ...
      @() pick(Q, @() gto2_gates01('gate0'), ...
                  @() gto2_gates01('gate0', 'nRef', 800, 'nTest', 30, 'nSeed', 1, 'hgrid', [1 0.25]));
  'S3',      'Gate 1: dimension excess, sigma sweep', ...
      @() pick(Q, @() gto2_gates01('gate1', 'nTest', 2500, 'seeds', 0:2, ...
                          'noiseGrid', [0.005 0.01 0.02 0.03 0.05 0.08 0.12 0.20 0.30]), ...
                  @() gto2_gates01('gate1', 'nRef', 800, 'nTest', 40, 'seeds', 0, ...
                          'noiseGrid', [0.01 0.30], 'dgrid', 2:3, 'nMC', 1e3));
  'S4',      'Gate 2: tangent rotation', ...
      @() pick(Q, @() gto2_paramsweep('rotation', 'nTest', 1000, 'seeds', 0:4), ...
                  @() gto2_paramsweep('rotation', 'nTest', 10, 'seeds', 0));
  'S5',      'Covariance route: curvature gradient', ...
      @() pick(Q, @() gto2_covroute('gradient'), ...
                  @() gto2_covroute('gradient', 'nTrial', 4, 'nRef', 1500, 'kgrid', [20 40]));
  'S6',      'Covariance route: boundary', ...
      @() pick(Q, @() gto2_covroute('boundary', 'kgrid', [20 40 80 160], 'nProbe', 200, 'nSeed', 3), ...
                  @() gto2_covroute('boundary', 'kgrid', 40, 'nProbe', 3, 'nSeed', 1, 'nRef', 2000, ...
                                    'fracs', [0 0.9]));
  'S7',      'Gate 3: Props 2-3 on quadric patches', ...
      @() pick(Q, @() gto2_curvlaw(), ...
                  @() gto2_curvlaw('nTrial', 4, 'nRef', 1000, 'kgrid', [20 40]));
  'S8',      'Gates 3-4: torus, convergence in nRef', ...
      @() pick(Q, @() torus_series([5000 10000 20000 40000]), ...
                  @() torus_series(3000, 'nTest', 400));
  'S9',      'Gate 3: anchor bias vs free parameters', ...
      @() pick(Q, @() gto2_paramsweep('anchor'), ...
                  @() gto2_paramsweep('anchor', 'nTrial', 3));
  'S10',     'Anisotropy vs the anchor (k = 40)', ...
      @() pick(Q, @() gto2_aniso('anchor'), ...
                  @() gto2_aniso('anchor', 'nTrialAnchor', 3));
  'S11',     'Anisotropy vs G,T,O (k = 100)', ...
      @() pick(Q, @() gto2_aniso('fractions'), ...
                  @() gto2_aniso('fractions', 'nRef', 600, 'nTest', 30, 'k', 40, 'kgrid', [0 400]));
  'S12',     'Anisotropy vs frame and tau^2 (k = 40)', ...
      @() gto2_aniso('frame');
  'S13',     'Anchor offset, pure on-manifold probes', ...
      @() pick(Q, @() gto2_anchorbias(), ...
                  @() gto2_anchorbias('nRef', 1500, 'nTest', 30, 'Rgrid', [3 20], 'kgrid', [20 40]));
  'S14',     'Gate 5: Clifford torus', ...
      @() pick(Q, @() gto2_manifolds('clifford'), ...
                  @() gto2_manifolds('clifford', 'nRef', 1500, 'nTest', 30, 'kgrid', [20 40]));
  'S15',     'Gate 5: circle x ellipse', ...
      @() pick(Q, @() gto2_manifolds('ellipse'), ...
                  @() gto2_manifolds('ellipse', 'nRefGrid', 3000, 'nTest', 100));
  'S16',     'Eq. 24 vs exact form vs anchor model', ...
      @() pick(Q, @() gto2_paramsweep('cor3'), ...
                  @() gto2_paramsweep('cor3', 'nTest', 30, 'seeds', 0, 'cells', [2000 40]));
  'S17',     'Off-manifold recovery, full grid', ...
      @() pick(Q, @() gto2_offmanifold(), ...
                  @() gto2_offmanifold('nRef', 1000, 'nTest', 30, 'kgrid', 40, 'nugrid', 0.3, ...
                                       'Ogrid', [0 0.5 1]));
  'S17fine', 'Off-manifold recovery, fine O* grid in the ill-conditioned block', ...
      @() pick(Q, @() gto2_offmanifold('signs', -1, 'kgrid', 160, 'nugrid', 0.3, 'Ogrid', 0:0.05:1), ...
                  @() gto2_offmanifold('signs', -1, 'kgrid', 40, 'nugrid', 0.3, 'Ogrid', [0 0.5], ...
                                       'nRef', 1000, 'nTest', 30));
};

R = struct();
R.meta = struct('stamp', stamp, 'quick', Q, 'logfile', logfile, 'version', version);
R.status = cell(0, 4);           % label, seconds, 'ok'/'FAILED', message
fprintf('\n%s\n  GTO2 RUN ALL  (%s)%s\n%s\n', repmat('#',1,78), stamp, ...
    ternary(Q, '  QUICK MODE', ''), repmat('#',1,78));

tAll = tic;
for i = 1:size(steps,1)
    label = steps{i,1};  desc = steps{i,2};  fh = steps{i,3};
    if ~isempty(o.only) && ~any(strcmp(label, o.only)), continue, end
    fprintf('\n%s\n  [%s] %s\n%s\n', repmat('=',1,78), label, desc, repmat('=',1,78));
    t = tic;
    try
        R.(label) = fh();
        R.status(end+1,:) = {label, toc(t), 'ok', ''};
    catch ME
        if exist('OCTAVE_VERSION', 'builtin')
            msg = ME.message;
            if ~isempty(ME.stack)
                msg = sprintf('%s (%s line %d)', msg, ME.stack(1).name, ME.stack(1).line);
            end
        else
            msg = getReport(ME, 'extended', 'hyperlinks', 'off');
        end
        fprintf(2, '\n  !! [%s] FAILED:\n%s\n', label, msg);
        R.status(end+1,:) = {label, toc(t), 'FAILED', msg};
    end
    save(matfile, 'R', '-v7');
end

fprintf('\n%s\n  SUMMARY (%.1f min total)\n%s\n', repmat('#',1,78), toc(tAll)/60, repmat('#',1,78));
for i = 1:size(R.status,1)
    fprintf('  %-8s %-7s %8.1f s\n', R.status{i,1}, R.status{i,3}, R.status{i,2});
end
nf = sum(strcmp(R.status(:,3), 'FAILED'));
fprintf('  %d step(s), %d failed.  Log: %s\n  Results: %s\n\n', size(R.status,1), nf, logfile, matfile);
save(matfile, 'R', '-v7');
end

% ---------------------------------------------------------------------
function out = pick(quick, fullFn, quickFn)
if quick, out = quickFn(); else, out = fullFn(); end
end

% ---------------------------------------------------------------------
function out = torus_series(nRefGrid, varargin)
out = struct('nRef', {}, 'torus', {});
for n = nRefGrid(:).'
    r = gto2_manifolds('torus', 'nRef', n, varargin{:});
    out(end+1).nRef = n;          %#ok<AGROW>
    out(end).torus  = r.torus;
end
% compact convergence summary: slope per bin at each nRef
fprintf('\n  Torus convergence: LS slope per theta bin (rows = bins, columns = nRef)\n');
fprintf('  %13s', 'theta bin');  fprintf(' %9d', nRefGrid);  fprintf('\n');
e = linspace(0, 2*pi, 9);
for b = 1:8
    fprintf('  %5.2f-%5.2f  ', e(b), e(b+1));
    for j = 1:numel(out)
        rb = out(j).torus.rows;
        idx = find(abs(rb(:,1) - e(b)) < 1e-9, 1);
        if isempty(idx), fprintf(' %9s', '--'); else, fprintf(' %9.3f', rb(idx,6)); end
    end
    fprintf('\n');
end
fprintf('\n');
end

% ---------------------------------------------------------------------
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
