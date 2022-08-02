function [allISIdist] = GSASmodel2(GSASparms,logtbins,numcells,numAS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%Note: vectr form uses logCV... for fitting
%logtbins should be base e

if ~isstruct(GSASparms)
    GSASparms = convertGSASparms(GSASparms,numcells,numAS);
end

%Ground state
allISIdist = LogGamma2(GSASparms.GSlogrates,GSASparms.GSCVs,GSASparms.GSweights,logtbins');
%%
%Add in each AS mode
%ASISI = zeros([size(allISIdist),length(GSASparms.ASlogrates)]);
for aa = 1:length(GSASparms.ASlogrates)
    allISIdist = allISIdist + LogGamma2(GSASparms.ASlogrates(aa),GSASparms.ASCVs(aa),GSASparms.ASweights(:,aa)',logtbins');
end
%%
%allISIdist = sum(ASISI,3)+GSISI;

end

