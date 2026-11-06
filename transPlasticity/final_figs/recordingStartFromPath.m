function d = recordingStartFromPath(f)
%RECORDINGSTARTFROMPATH  datenum of recording start from a session folder name
%   such as CP12_250604_075757 (YYMMDD_HHMMSS). NaN if there is no HHMMSS token.
d = NaN;
tok = regexp(f, '_(\d{6})_(\d{6})(?:\\|/|$)', 'tokens');
if isempty(tok), return; end
s = tok{end}{1}; h = tok{end}{2};
yy = str2double(s(1:2)); mm = str2double(s(3:4)); dd = str2double(s(5:6));
hh = str2double(h(1:2)); mi = str2double(h(3:4)); ss = str2double(h(5:6));
if mm < 1 || mm > 12 || dd < 1 || dd > 31 || hh > 23 || mi > 59 || ss > 59, return; end
d = datenum(2000 + yy, mm, dd, hh, mi, ss);
end
