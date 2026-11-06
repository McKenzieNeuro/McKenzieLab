function ref = firstASVReference(fils, ani, cnd, dnm, inf_t, recStart, expCond)
%FIRSTASVREFERENCE  Per animal, the first ASV infusion (session, date, clock, time).
%
%   ref = firstASVReference(fils, ani, cnd, dnm, inf_t, recStart, expCond)
%
%   fils     cellstr of session files        ani, cnd  cellstr animal / condition
%   dnm      datenum of the session date     inf_t     first Infusion_start (s), NaN if none
%   recStart datenum of recording start (NaN when the folder name has no HHMMSS)
%   expCond  condition string of exposure sessions, e.g. 'ASV'
%
%   ref.animal  cellstr;  ref.item(k) fields: session, date, recStart, infusionSec, clock.
%   The earliest exposure session with a marked infusion wins (by clock time when
%   known, otherwise by date). See TS2FROMREFERENCE.
ua = unique(ani(:));
blank = struct('session','','date',NaN,'recStart',NaN,'infusionSec',NaN,'clock',false);
ref = struct('animal',{ua}, 'item',repmat(blank, numel(ua), 1));
for u = 1:numel(ua)
    m = find(strcmp(ani, ua{u}) & strcmp(cnd, expCond) & isfinite(inf_t));
    if isempty(m), continue; end
    key = dnm(m) + 1e-6*(1:numel(m))';
    key(isfinite(recStart(m))) = recStart(m(isfinite(recStart(m))));
    [~, o] = min(key); j = m(o);
    ref.item(u) = struct('session',fileparts(fils{j}), 'date',dnm(j), 'recStart',recStart(j), ...
        'infusionSec',inf_t(j), 'clock',isfinite(recStart(j)));
end
end
