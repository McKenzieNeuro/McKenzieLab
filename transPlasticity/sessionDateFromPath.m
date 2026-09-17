function dn = sessionDateFromPath(f)
%SESSIONDATEFROMPATH  Datenum of a session from its folder name.
%
%   dn = sessionDateFromPath(f)
%
%   Expects a session folder like  CP12_250604  or  CP12_250604_075757
%   (animal, then YYMMDD, optionally then HHMMSS) anywhere in the path, and
%   returns MATLAB datenum. NaN if no 6-digit date-like token is found.
%
%   The 6-digit group is read as YYMMDD with a 2000s century. Sanity-checked
%   on month and day: a token that parses to month 00 or 13+, or day 00 or
%   32+, is rejected rather than silently producing a nonsense date, because
%   a mis-parsed date here would silently corrupt the pre/post ordering that
%   the whole analysis rests on.
%
%   Uses the LAST matching token in the path, so a parent folder that happens
%   to contain digits does not win over the session folder itself.

dn = NaN;
tok = regexp(f, '_(\d{6})(?:_|\\|/|$)', 'tokens');
if isempty(tok), return; end

for k = numel(tok):-1:1          % last match first
    s = tok{k}{1};
    yy = str2double(s(1:2));
    mm = str2double(s(3:4));
    dd = str2double(s(5:6));
    if mm >= 1 && mm <= 12 && dd >= 1 && dd <= 31
        dn = datenum(2000 + yy, mm, dd);
        return
    end
end
end
