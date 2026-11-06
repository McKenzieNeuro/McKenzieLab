function m = spikeDistanceBinMean(E, y1, y2, g)
%SPIKEDISTANCEBINMEAN  Exact mean of a piecewise-linear profile over bins.
%
%   m = spikeDistanceBinMean(E, y1, y2, g)
%
%   E, y1, y2 from SPIKEDISTANCEPROFILE; g = bin edges (sorted, inside
%   [E(1) E(end)]). m(b) = (1/(g(b+1)-g(b))) * integral of S(t) over the bin,
%   computed exactly from the linear pieces (no resampling).
E = E(:); y1 = y1(:); y2 = y2(:); g = g(:);
F = [0; cumsum(0.5 * (y1 + y2) .* diff(E))];         % integral up to each E
j = discretize(g, E);                                 % interval containing each edge
h = E(j+1) - E(j);
dx = g - E(j);
Fg = F(j) + y1(j) .* dx + (y2(j) - y1(j)) .* dx.^2 ./ (2 * h);
m = diff(Fg) ./ diff(g);
end
