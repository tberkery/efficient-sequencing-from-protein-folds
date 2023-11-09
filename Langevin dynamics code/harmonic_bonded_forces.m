function forces = harmonic_bonded_forces(pos1, pos2, k, le)

xyz_bond = pos1-pos2;
r_bond = sqrt(xyz_bond(1)^2 + xyz_bond(2)^2 + xyz_bond(3)^2);

% Harmonic bonded potential equation
bracket = -k * (r_bond - le) / r_bond;
forces = xyz_bond * bracket;

end