function [loc, how] = matchCellMetrics(cm, spikesIdx, cluID, nSpikeUnits)
%MATCHCELLMETRICS  Row of cell_metrics for each unit.
%
%   [loc, how] = matchCellMetrics(cm, spikesIdx, cluID, nSpikeUnits)
%
%   spikesIdx    index of each unit into spikes.times
%   cluID        that unit's cluster ID (spikes.UID in these files holds the
%                sorter's cluster IDs, e.g. 0 5 8 10 30), [] if unknown
%   nSpikeUnits  numel(spikes.times), [] if unknown
%   loc          index into cell_metrics.* (0 = no match);  how = 'index'|'cluID'|'none'
%
%   cell_metrics.UID is 1..N (CellExplorer renumbers), so it must NOT be
%   matched against spikes.UID. CellExplorer keeps the spikes order, so the
%   spikes index is used when the unit counts agree, and checked against
%   cell_metrics.cluID when both are available; on disagreement cluID wins.
spikesIdx = spikesIdx(:); cluID = cluID(:);
nC = numel(cm.putativeCellType);
loc = zeros(numel(spikesIdx), 1); how = 'none';
hasClu = isfield(cm, 'cluID') && numel(cm.cluID) == nC && ~isempty(cluID) && numel(cluID) == numel(spikesIdx);
sameN = isempty(nSpikeUnits) || nSpikeUnits == nC;
if sameN && all(spikesIdx >= 1 & spikesIdx <= nC)
    loc = spikesIdx; how = 'index';
    if hasClu
        cc = cm.cluID(:);
        bad = cc(loc) ~= cluID;
        if any(bad)
            [tf, l] = ismember(cluID, cc);
            if all(tf)          % the IDs are cluster IDs and the order differs: trust cluID
                warning('matchCM:clu', '%d of %d units: index order disagrees with cell_metrics.cluID; matching on cluID', ...
                    sum(bad), numel(bad));
                loc = l; how = 'cluID';
            end                 % otherwise the given IDs are not cluster IDs: keep index order
        end
    end
elseif hasClu
    [tf, l] = ismember(cluID, cm.cluID(:)); loc(tf) = l(tf); how = 'cluID';
end
end
