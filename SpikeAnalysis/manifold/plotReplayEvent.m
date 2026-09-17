function plotReplayEvent(out, eventID)
% plotReplayEvent  Visualize one decoded replay event from bayesianReplayAnalysis.
%
%   plotReplayEvent(out, eventID)
%
%   out     : output struct from bayesianReplayAnalysis.
%   eventID : row index into out.results / out.posteriors.

P = out.posteriors{eventID};
if isempty(P)
    warning('Event %d was excluded or posteriors were not saved.', eventID);
    return;
end
bc  = out.binCenters;
tau = out.opts.tau;
r   = out.results(eventID, :);
[nPos, nT] = size(P);
t = (0:nT-1) * tau * 1000;          % ms

imagesc(t, bc, P); axis xy; hold on;
colormap(hot); cb = colorbar; cb.Label.String = 'P(position | spikes)';
xlabel('Time in event (ms)'); ylabel('Linearized position (cm)');

% overlay best-fit line if available
S = scoreTrajectory(P, bc, tau, struct('computeLine', true, ...
        'lineBand', out.opts.lineBand, 'nSlope', out.opts.nSlope));
if isfield(S,'lineParams') && all(isfinite(S.lineParams))
    a = S.lineParams(1); b = S.lineParams(2);
    j = 0:nT-1;
    yk = bc(min(max(round(a*j + b),1),nPos));
    plot(j*tau*1000, yk, 'c-', 'LineWidth', 2);
end

title(sprintf(['Event %d | wcorr=%.2f  line=%.2f  slope=%.0f cm/s  ' ...
    'p=%.3f (%s)  cells=%d'], r.eventID, r.weightedCorr, r.lineScore, ...
    r.lineSlope_cmps, r.pValue, r.scoreType{1}, r.nActiveCells));
hold off;
end
