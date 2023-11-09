% Calculate force in x y z row vector that particle j exerts on j
% Requires file force_calc

function Fij = lennard_force_noMIC(ri, rj, epsilon, sigma)

r = ri-rj;
r_sq = r(1)^2 + r(2)^2 + r(3)^2;

% Calculate force ij for x,y,and z
inner_bracket = (48*epsilon/r_sq) * ( (sigma^2/r_sq)^6 - 1/2 * (sigma^2/r_sq)^3 );
Fij = zeros(1,3); % Initialize Fij
Fij(1) = inner_bracket*r(1);
Fij(2) = inner_bracket*r(2);
Fij(3) = inner_bracket*r(3);

end