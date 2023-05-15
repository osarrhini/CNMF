function ptac = plasma_from_blood_activity(btac,t)
%%blood_to_plasma Plasma activty from whole blood activity as a function of
% time.
% Inputs:
%   btac: Whole blood concentration of 18F-FDG
%   t: Time must be in seconds.
% Output:
%   ptac: Plasma concentration of 18F-FDG
% See:
% J Nucl Med 2007; 48:837–845
% DOI: 10.2967/jnumed.106.038182
%
% To call: ptac = plasma_from_blood_activity(btac,t);
% @Author: Otman Sarrhini, Sherbrooke University
% **This code comes with no guarantee or warranty of any kind.**

ptac=btac(:).*(0.386*exp(-0.191*t(:)/60)+1.165);