function [csd, CSDelecinds] = trad_csd(xcpots, elec_sep_mm, xc_cond, FORCE)
% TRAD_CSD [csd, CSDelecinds] = trad_csd(xcpots, elec_sep_mm, xc_cond*, FORCE*)
% (* = optional argument)
%
% Traditional estimation of the CSD using linear electrode recordings.
% FORCE flag prevents function from treating larger dimension as time
% dimension. If xc_cond is set to 1, csd is returned in units of mV/mm^2.

if (nargin < 3) || isempty(xc_cond) || (xc_cond <= 0)
	xc_cond = 1/3333; % S/mm = 3e-3 S/cm
end
if (nargin < 4)
    FORCE = 0;
end

if (size(xcpots,2) > size(xcpots,1)) && ~FORCE
    xcpots = xcpots';
end
N = size(xcpots,2);

% There also must be at least three electrodes
if (N < 3)
    error('Voltages from at least three electrodes needed')
end

% Create 2nd order derivative matrix; for N electrodes, this is size (N-2) x N.
% The square of the electrode separation is in the denominator.
D = zeros(N-2,N);
D(1,1:3) = [1 -2 1]/elec_sep_mm^2; % mm^-2
for i=2:N-2
    D(i,:) = circshift(D(i-1,:),[0,1]);
end

% CSD = -(conductivity)*(2nd spatial derivative of voltage along electrode line)
CSDelecinds = 2:(size(xcpots,2)-1);
csd = (-xc_cond*D*xcpots')'; % (ohm-mm)^-1*(mV/mm^2) = mA/mm^3
