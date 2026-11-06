function [ts2, mode, dayOff, rf] = ts2FromReference(ts1, animal, date_, rs, ref)
%TS2FROMREFERENCE  Seconds from the animal's first ASV infusion for times ts1
%   (s from the start of this recording).
%
%   [ts2, mode, dayOff, rf] = ts2FromReference(ts1, animal, date_, rs, ref)
%
%   date_  datenum of this session's date;  rs  datenum of this recording's start
%   (NaN if the folder name has no HHMMSS);  ref  from FIRSTASVREFERENCE.
%
%   mode 'clock'  exact: both recordings have a start clock time.
%        'date'   only dates are known: ts2 = days*86400 + (ts1 - first-ASV
%                 infusion time), which assumes this recording started at the
%                 same time of day as the first-ASV one (error up to +-24 h).
%        'none'   the animal has no ASV session; ts2 = NaN.
%   dayOff is the exact integer day offset from the first-ASV session.
rf = struct('session','','date',NaN,'recStart',NaN,'infusionSec',NaN,'clock',false);
k = find(strcmp(ref.animal, animal), 1);
if ~isempty(k), rf = ref.item(k); end
ts2 = nan(size(ts1)); mode = 'none'; dayOff = NaN;
if ~isempty(rf.session) && isfinite(rf.infusionSec)
    dayOff = round(date_ - rf.date);
    if isfinite(rs) && rf.clock
        ts2 = ((rs + ts1/86400) - (rf.recStart + rf.infusionSec/86400)) * 86400;
        mode = 'clock';
    else
        ts2 = dayOff*86400 + (ts1 - rf.infusionSec);
        mode = 'date';
    end
end
end
