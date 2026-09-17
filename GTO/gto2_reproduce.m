function out = gto2_reproduce(which)
% GTO2_REPRODUCE  Regenerate every table in the supplement, in supplement order.
%
%   gto2_reproduce            % everything (slow: the ellipse and gate-0 sweeps dominate)
%   gto2_reproduce('fast')    % everything except the four heaviest sweeps
%   gto2_reproduce('checks')  % the invariant checks C1-C5 only
%   gto2_reproduce('S7')      % a single table by its supplement number
%
% Each block prints a banner naming the supplement table it reproduces, then calls
% the same entry point the supplement's caption names. Nothing here computes
% anything itself: this file exists so that the mapping from table to code is
% executable rather than a claim in a caption.
%
% TABLE MAP (supplement numbering)
%   S1  exact identity, Prop 3          gto2_checks('c5')
%   S2  Corollary 3 residual            gto2_paramsweep('cor3')
%   S3  gate 0, sheet proximity         gto2_gates01('gate0')
%   S4  gate 1, dimension excess        gto2_gates01('gate1')
%   S5  gate 2, frame rotation          gto2_paramsweep('rotation')
%   S6  gate 2, curvature gradient      gto2_covroute('gradient')
%   S7  gate 2, boundary                gto2_covroute('boundary')
%   S8  gate 3, anchor displacement law gto2_curvlaw
%   S9  gate 3, free parameters         gto2_paramsweep('anchor')
%   S10 gate 3/4, torus of revolution   gto2_manifolds('torus')
%   S11 gate 4, anisotropy vs fractions gto2_aniso('fractions')
%   S12 gate 3, anisotropy vs anchor    gto2_aniso('anchor')
%   S13 gate 2/6, anisotropy vs frame   gto2_aniso('frame')
%   S14 gate 5, booking at |c| = 1      gto2_anchorbias
%   S15 gate 5, booking at |c| < 1      gto2_manifolds('clifford')
%   S16 gate 5, varying second form     gto2_manifolds('ellipse')
%   S17 gate 6, recovery of O-hat       gto2_offmanifold
%   S18 cascade summary                 derived from the above; not generated
%
% ON REPRODUCING THE PUBLISHED NUMBERS. The tables in the supplement were produced
% by a reference implementation, not by this suite. Expect agreement in every
% scaling, sign and ordering, and expect the third decimal to move: the random
% streams differ. Where a published value is quoted to more figures than the
% seed-to-seed spread supports, the spread is the thing to trust. Any DISAGREEMENT
% IN SIGN OR SCALING is a real discrepancy and should be treated as one.
%
% RETURNS out.<field> for each block run, keyed by table number.

if nargin < 1, which = 'all'; end
which = lower(which);
out = struct();

blocks = { ...
 'S1',  'exact identity (Prop 3)',            @() gto2_checks('c5'),            false; ...
 'S2',  'Corollary 3 residual',               @() gto2_paramsweep('cor3'),      false; ...
 'S3',  'gate 0: sheet proximity',            @() gto2_gates01('gate0'),        true ; ...
 'S4',  'gate 1: dimension excess',           @() gto2_gates01('gate1'),        false; ...
 'S5',  'gate 2: frame rotation',             @() gto2_paramsweep('rotation'),  false; ...
 'S6',  'gate 2: curvature gradient',         @() gto2_covroute('gradient'),    false; ...
 'S7',  'gate 2: boundary',                   @() gto2_covroute('boundary'),    false; ...
 'S8',  'gate 3: anchor displacement law',    @() gto2_curvlaw,                 false; ...
 'S9',  'gate 3: free-parameter dependence',  @() gto2_paramsweep('anchor'),    true ; ...
 'S10', 'gates 3 and 4: torus of revolution', @() gto2_manifolds('torus'),      false; ...
 'S11', 'gate 4: anisotropy, fractions',      @() gto2_aniso('fractions'),      false; ...
 'S12', 'gate 3: anisotropy, anchor',         @() gto2_aniso('anchor'),         true ; ...
 'S13', 'gates 2 and 6: anisotropy, frame',   @() gto2_aniso('frame'),          false; ...
 'S14', 'gate 5: booking at |c| = 1',         @() gto2_anchorbias,              false; ...
 'S15', 'gate 5: booking at |c| < 1',         @() gto2_manifolds('clifford'),   false; ...
 'S16', 'gate 5: varying second form',        @() gto2_manifolds('ellipse'),    true ; ...
 'S17', 'gate 6: recovery of O-hat',          @() gto2_offmanifold,             false  };

switch which
    case 'checks'
        fprintf('\n%s\n  INVARIANT CHECKS C1-C5\n%s\n', repmat('=',1,78), repmat('=',1,78));
        out.checks = gto2_checks('all');
        return
    case {'all','fast'}
        sel = 1:size(blocks,1);
        if strcmp(which,'fast')
            sel = sel(~[blocks{:,4}]);
            fprintf('\n  [fast] skipping the four heaviest sweeps: S3, S9, S12, S16\n');
        end
    otherwise
        sel = find(strcmpi(blocks(:,1), which));
        if isempty(sel)
            error('gto2_reproduce:table', ...
                'Unknown selector "%s". Use all, fast, checks, or a table id S1-S17.', which);
        end
end

t0 = tic;
for i = sel(:).'
    fprintf('\n%s\n  TABLE %s  --  %s\n  call: %s\n%s\n', ...
        repmat('=',1,78), blocks{i,1}, blocks{i,2}, func2str(blocks{i,3}), repmat('=',1,78));
    try
        out.(blocks{i,1}) = blocks{i,3}();
    catch err
        fprintf(2, '  *** %s FAILED: %s\n', blocks{i,1}, err.message);
        out.(blocks{i,1}) = err;
    end
end
fprintf('\n%s\n  done in %.1f s\n', repmat('-',1,78), toc(t0));

failed = {};
f = fieldnames(out);
for i = 1:numel(f)
    if isa(out.(f{i}), 'MException'), failed{end+1} = f{i}; end %#ok<AGROW>
end
if isempty(failed)
    fprintf('  all %d blocks returned\n\n', numel(sel));
else
    fprintf(2, '  %d block(s) failed: %s\n\n', numel(failed), strjoin(failed, ', '));
end
end
