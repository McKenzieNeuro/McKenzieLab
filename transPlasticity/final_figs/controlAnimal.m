function a = controlAnimal(f)
%CONTROLANIMAL  Animal ID of a control (Huszar) session from its file path:
%   the session folder name up to its first _YYMMDD token,
%   e.g. R:\HuszarR\e13_26m1_211019\x.mat -> e13_26m1.
[~, s] = fileparts(fileparts(f));
a = regexprep(s, '_\d{6}.*$', '');
end
