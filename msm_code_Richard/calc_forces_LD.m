%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Melissa Mai
% March 2017
%
% Function calc_force
% Calculates the pairwise forces between all particles and
% the potential energy of the system according to the Lennard-Jones
% Potential for MD simulations.
%
% Force and energy calculations are adapted from Algorithms 5 and 6 from
% Chapter 4 of Frenkel & Smit's Understanding Molecular Simulations.
%
% INPUTS
% coords: coordinates for the particles in the system
% box: box length
% sigma: distance where interparticle potential is zero
% eps: energy scaling factor in LJ potential
%
% OUTPUTS
% forces: matrix of the total force on each particle in the x, y, and z
% directions
% enpot: Total Lennard-Jones potential energy of the system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [forces, PE] = calc_forces_LD(hpChain, coords, velocities, dt, temp, mass, friction, k, le)
    [num_particles, num_elements] = size(coords);
    [forces_config, PE_config] = calc_forces_configuration(hpChain, coords, k, le);
    forces_random = randn(num_particles, num_elements) * sqrt(2*temp*mass*friction/dt);
    forces_dissipative = friction*mass.*velocities;
    forces = forces_config + forces_random - forces_dissipative;
    PE = PE_config;
end