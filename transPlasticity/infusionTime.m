function t = infusionTime(spikeFile, whichOne)
%INFUSIONTIME  Infusion start time (s) for a session, from its .evt.sti file.
%
%   t = infusionTime(spikeFile)            -> first Infusion_start, or NaN
%   t = infusionTime(spikeFile, 'all')     -> all Infusion_start times
%
%   spikeFile is the *.spikes.cellinfo.mat path; the event file is looked up
%   in the same folder as amplifier_analogin_auxiliary_int16.evt.sti.
%
%   Returns NaN (with no error) when there is no event file or no
%   Infusion_start marker -- SAL/Baseline sessions may legitimately have
%   neither, and callers need to distinguish "no infusion here" from a hard
%   failure. A session that HAS an event file but no matching description is
%   warned about, since that more often means the description string differs
%   from what is expected than that the infusion truly did not happen.
%
%   Units are whatever LoadEvents returns. Check against your session
%   duration before trusting it: if the value looks like milliseconds or
%   samples rather than seconds, everything downstream that compares it to
%   spike times will silently mis-split the session.
%
%   See also SYNCBATCH, PREPOSTASV.

if nargin < 2, whichOne = 'first'; end
t = NaN;

d = fileparts(spikeFile);
evFile = fullfile(d, 'amplifier_analogin_auxiliary_int16.evt.sti');
if ~exist(evFile, 'file')
    c = dir(fullfile(d, '*.evt.sti'));
    if isempty(c), return; end
    evFile = fullfile(c(1).folder, c(1).name);
end

try
    ev = LoadEvents(evFile);
catch ME
    warning('infusionTime:load','%s: LoadEvents failed (%s)', evFile, ME.message);
    return
end

if ~isfield(ev,'time') || ~isfield(ev,'description'), return; end
hit = contains(ev.description, 'Infusion_start');
if ~any(hit)
    warning('infusionTime:noMarker', ...
        ['%s: event file present but no ''Infusion_start'' description ' ...
         '(found: %s). Check the marker string.'], evFile, ...
        strjoin(unique(ev.description(:))', ', '));
    return
end

tt = ev.time(hit);
if strcmpi(whichOne, 'all')
    t = tt(:);
else
    t = min(tt);
end
end
