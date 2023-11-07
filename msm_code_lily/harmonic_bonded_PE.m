% Harmonic bonded potential. Calculate potential energy of a stretching
% covalent bond. pos1 and pos2 are the 3D positions of two bonded particles
function potential = harmonic_bonded_PE(pos1, pos2, k, le)

% Calculate distance between the particles
xyz_bond = pos1-pos2;
r_bond = sqrt(xyz_bond(1)^2 + xyz_bond(2)^2 + xyz_bond(3)^2);

% Harmonic bonded potential equation
potential = 0.5 * k * (r_bond - le)^2;

end