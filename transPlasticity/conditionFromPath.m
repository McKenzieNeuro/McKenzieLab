function cond = conditionFromPath(f, rootName)
%CONDITIONFROMPATH  Condition = the path component two levels below rootName.
%
%   cond = conditionFromPath(f, rootName)
%
%   R:\TransPlasticity\CP19\ASV24\CP19_250604\...  + 'TransPlasticity' -> 'ASV24'
%
%   Layout assumed:  <root>\<animal>\<condition>\<session>\<file>
%   i.e. condition sits immediately below the animal folder that
%   animalFromPath returns. Returns '' if the path is shorter than that or
%   the root isn't found -- check for empties rather than assuming every
%   path parses, since a session stored one level shallower or deeper will
%   silently return the wrong folder name otherwise.
%
%   See also ANIMALFROMPATH.

parts = regexp(f, '[\\/]', 'split');
parts = parts(~cellfun(@isempty, parts));
k = find(strcmpi(parts, rootName), 1, 'first');
if isempty(k) || k+2 > numel(parts)
    cond = '';
    return
end
cond = parts{k+2};
end