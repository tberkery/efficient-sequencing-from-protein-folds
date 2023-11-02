%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Melissa Mai
% March 2017
%
% function initialize_velocities
% Initializes velocities for an MD simulation. This function generates
% normally distributed random numbers for each particle in the x, y, and z
% directions. The average velocity is subtracted from each term so that the
% total momentum of the system is 0. The velocities are then scaled to 
% achieve the desired temperature.
% 
% Units are reduced.
% i.e., temp = Temp * kb / eps
%
% INPUTS
% npart: number of particles
% temp: temperature
%
% OUTPUTS
% vels: matrix of velocities for each particle in the x, y, and z
% directions
% actualtemp: actual temperature (reduced units) given generated initial
% velocities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vels] = initialize_velocities(npart, temp)
    % Generate normally distributed velocities for all particles in all
    % directions
    vels = randn(npart, 3) * sqrt(temp);
    avg_vels = sum(vels) / npart;
    
    % Set momentum to zero
    vels = vels - avg_vels;
end