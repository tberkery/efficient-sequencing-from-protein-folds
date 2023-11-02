%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Michaelis-Menten model of product from input parameters, substrate. 
% 
% Parameters:
%   beta
%       beta(1) - V_max, the max velocity of the reaction
%       beta(2) - K_m, the Michaelis constant
%   x - vector of substrate concentrations
% 
% Return:
%   y - vector of reaction velocities (modeled by Michaelis-Menten)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[y]=MM(beta, x)
    V_max = beta(1);
    K_m = beta(2);
    y = (V_max * x) ./ (K_m + x);
end