function id = animalFromPath(f, rootName)
%ANIMALFROMPATH  Animal ID = the path component immediately after rootName.
%
%   id = animalFromPath(f, rootName)
%
%   R:\TransPlasticity\CP12\ASV\CP12_250604\...   + 'TransPlasticity' -> 'CP12'
%   R:\WSun\ePhys\OptoRAM\OptoRAM7\OptoRAM7_...   + 'OptoRAM'         -> 'OptoRAM7'
%
%   Matching is exact and case-insensitive on the folder name, and takes the
%   FIRST match, so a root of 'OptoRAM' resolves against the 'OptoRAM' folder
%   and returns the 'OptoRAM7' below it rather than matching 'OptoRAM7'
%   itself. Returns '' if the root isn't found or nothing follows it -- check
%   for empties rather than assuming every path parses.

parts = regexp(f, '[\\/]', 'split');
parts = parts(~cellfun(@isempty, parts));
k = find(strcmpi(parts, rootName), 1, 'first');
if isempty(k) || k >= numel(parts)
    id = '';
    return
end
id = parts{k+1};
end
