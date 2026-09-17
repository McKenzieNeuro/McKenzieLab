%% demo_bayesianReplay.m
% Bayesian replay analysis over ripple-defined candidate events.
% Files: computeRateMaps1D.m, bayesianDecode.m, scoreTrajectory.m,
%        bayesianReplayAnalysis.m, plotReplayEvent.m
%
% --------------------------------------------------------------------------
% A) REAL DATA (your Buzcode session) — template, edit paths/fields and run.
% --------------------------------------------------------------------------
% basename = 'RO1_240106';
% load([basename '.spikes.cellinfo.mat']);     % -> spikes  (spikes.times{:})
% load([basename '.Behavior.mat']);            % -> behavior (.timestamps, .position.lin)
% load([basename '.ripples.events.mat']);      % -> ripples  (.timestamps Nx2)
% load([basename '.session.mat']);             % -> session (.epochs for run window)
% % If you have it: load([basename '.cell_type.mat']); -> cell_type (cellstr)
%
% % Restrict place-field construction to the RUN/track epoch:
% sessionID = 2;
% runEpoch  = [session.epochs{sessionID}.startTime, session.epochs{sessionID}.stopTime];
%
% % Use pyramidal cells only (replay decoding uses place cells):
% if exist('cell_type','var')
%     cellIdx = find(strcmpi(cell_type, 'Pyramidal Cell'));
% else
%     cellIdx = 1:numel(spikes.times);
% end
%
% opts = struct( ...
%     'cellIdx',       cellIdx, ...
%     'runEpoch',      runEpoch, ...
%     'tau',           0.02, ...     % 20 ms decoding bins
%     'binSize',       2, ...        % cm  (make sure behavior.position.lin is in cm!)
%     'smoothSD',      4, ...        % cm
%     'speedThresh',   5, ...        % cm/s for "running"
%     'minActiveCells',5, ...        % event inclusion
%     'minBins',       5, ...
%     'minSpkBins',    5, ...
%     'peakWindow',    0.1, ...      % decode a fixed 100 ms window around each
%                                    % ripple peak (5 bins at tau=20 ms). Use []
%                                    % to decode the raw detected ripple bounds.
%     'scoreForP',     'wcorr', ...  % 'wcorr' (fast) or 'line' (radon)
%     'shuffleMethod', 'columnCycle', ... % or 'permute'
%     'nShuffles',     500);
%
% out = bayesianReplayAnalysis(spikes, behavior, ripples, opts);
% sig = out.results(out.results.included & out.results.pValue < 0.05, :);
% disp(sortrows(sig, 'pValue'));
% figure; plotReplayEvent(out, sig.eventID(1));     % look at the top event
%
% Notes:
%  * behavior.position.lin must be in the SAME units as binSize/speedThresh (cm).
%    The RippleTagging/Yang-2024 loader exposes it as behavior.position.lin.
%  * If your spike times come from a *.spike_time.mat cell (as in the UMAP demo),
%    pass that cell directly as the first argument instead of a spikes struct:
%        out = bayesianReplayAnalysis(spike_time, behavior, ripples, opts);
%  * Linear-track caveat: place fields are direction-dependent. For a true linear
%    track, build/decode the two travel directions separately (split runEpoch or
%    pass per-direction rateMaps via opts.rateMaps) and keep the better-scoring
%    direction per event.

%% --------------------------------------------------------------------------
% B) SELF-TEST on synthetic data (runnable now; verifies the whole pipeline).
% --------------------------------------------------------------------------
RUN_SELFTEST = true;
if RUN_SELFTEST
    rng(1);

    % ---- ground-truth 1D place fields ----
    binSize = 2; trackLen = 160;
    edges = 0:binSize:trackLen; xc = edges(1:end-1) + binSize/2; nPos = numel(xc);
    nCells = 60;
    centers = linspace(5, trackLen-5, nCells);
    Rtrue = 12 * exp(-(xc - centers(:)).^2 ./ (2*8^2)) + 0.05;   % nCells x nPos (Hz)

    % ---- synthetic running behavior (triangle sweeps of the track) ----
    dt = 0.025; Tdur = 240; posTime = (0:dt:Tdur).';
    sweep = 25;                                   % cm/s
    lin = trackLen/2 + (trackLen/2) * sawtooth_(2*pi*(sweep/(2*trackLen))*posTime, 0.5);
    lin = min(max(lin, 0), trackLen);
    behavior = struct('timestamps', posTime, 'position', struct('lin', lin));

    % ---- run spikes from the true fields (so fields can be rebuilt) ----
    spikeTimes = cell(nCells, 1);
    posbinAll = max(1, min(nPos, ceil(lin / binSize)));
    for c = 1:nCells
        lam = Rtrue(c, posbinAll).' * dt;          % expected count per sample
        k = poissSample(lam);
        idx = find(k > 0);
        st = [];
        for ii = idx.'
            st = [st; posTime(ii) + dt*rand(k(ii),1)]; %#ok<AGROW>
        end
        spikeTimes{c} = sort(st);
    end

    % ---- candidate events placed AFTER the run epoch ----
    tau = 0.02; nEvt = 40; evDur = 0.2; gap = 0.5; t0 = Tdur + 2;
    evt = zeros(nEvt, 2); isTraj = false(nEvt, 1);
    for e = 1:nEvt
        s = t0 + (e-1)*(evDur + gap);
        evt(e,:) = [s, s + evDur];
        isTraj(e) = rand < 0.5;                    % half real trajectories, half null
        nb = round(evDur / tau);
        if isTraj(e)
            p0 = rand*0.3*trackLen + 0.1*trackLen;
            v  = (400 + rand*500) * sign(rand-0.5); % cm/s replay speed
            gain = 6;
            for b = 1:nb
                pos = p0 + v*(b-1)*tau;
                pb  = max(1, min(nPos, round(pos/binSize)));
                k = poissSample(Rtrue(:,pb)*tau*gain);
                for c = find(k>0).'
                    spikeTimes{c} = [spikeTimes{c}; s + (b-1)*tau + tau*rand(k(c),1)];
                end
            end
        else
            base = mean(Rtrue(:))*tau*6;
            for b = 1:nb
                k = poissSample(base*ones(nCells,1));
                for c = find(k>0).'
                    spikeTimes{c} = [spikeTimes{c}; s + (b-1)*tau + tau*rand(k(c),1)];
                end
            end
        end
    end
    for c = 1:nCells, spikeTimes{c} = sort(spikeTimes{c}); end
    ripples = struct('timestamps', evt);

    % ---- run the pipeline (build fields from run epoch, decode events) ----
    opts = struct('runEpoch',[0 Tdur], 'tau',tau, 'binSize',binSize, ...
                  'smoothSD',4, 'speedThresh',5, 'scoreForP','wcorr', ...
                  'shuffleMethod','columnCycle', 'nShuffles',500, 'verbose',true);
    out = bayesianReplayAnalysis(spikeTimes, behavior, ripples, opts);

    % ---- field reconstruction sanity (size-robust, no toolbox) ----
    [~, bi] = max(out.rateMaps, [], 2); peakRebuilt = out.binCenters(bi);
    [~, bt] = max(Rtrue, [], 2);        peakTrue    = xc(bt);
    medErr = median(abs(peakRebuilt(:) - peakTrue(:)));
    fprintf('\nMedian place-field peak error (rebuilt vs true): %.1f cm\n', medErr);

    % ---- separation of planted trajectories vs null ----
    R = out.results; inc = R.included;
    fprintf('Median p  | trajectory events: %.3f   null events: %.3f\n', ...
        median(R.pValue(inc & isTraj)), median(R.pValue(inc & ~isTraj)));
    fprintf('Frac p<0.05 | trajectory: %.2f   null: %.2f\n', ...
        mean(R.pValue(inc & isTraj) < 0.05), mean(R.pValue(inc & ~isTraj) < 0.05));

    % ---- show one significant trajectory event ----
    cand = find(inc & isTraj & R.pValue < 0.05);
    if ~isempty(cand)
        figure; plotReplayEvent(out, cand(1));
    end
end

%% ---- local helpers used only by the self-test ----
function k = poissSample(L)
% Vectorized Knuth Poisson sampler (no Statistics Toolbox needed).
sz = size(L); L = L(:);
Lth = exp(-L);
p = ones(numel(L),1); k = zeros(numel(L),1);
notdone = true(numel(L),1);
while any(notdone)
    p(notdone) = p(notdone) .* rand(nnz(notdone),1);
    inc = notdone & (p > Lth);
    k(inc) = k(inc) + 1;
    notdone = inc;
end
k = reshape(k, sz);
end

function y = sawtooth_(t, width)
% Minimal triangle/sawtooth (avoids Signal Processing Toolbox).
% width=0.5 gives a symmetric triangle wave in [-1,1].
t = mod(t/(2*pi), 1);
y = zeros(size(t));
up = t < width;
y(up)  = -1 + 2*t(up)/width;
y(~up) =  1 - 2*(t(~up)-width)/(1-width);
end
