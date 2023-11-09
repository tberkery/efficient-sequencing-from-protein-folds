% Force dissipation on bead motion using the low friction Langevin form
% Force in one direction of x y or z. zeta is friction coefficient.

function Fx = dissipative_force(vx, m, zeta)
R_debug = 1;
Fx = R_debug * (-zeta) * m * vx;
end