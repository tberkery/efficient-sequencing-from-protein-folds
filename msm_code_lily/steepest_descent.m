% Steepest Descent for homework 6
% init_coords is the initial coordinates of all particles in a Nx3 array
% hp is the hydrophobicity of N particles in a Nx1 array
% sigma is the threshold energy difference to stop the minimization
% init_lambda is the initial step size 

function [min_coords, minE] = steepest_descent(init_coords,hp,sigma,stepsize,dec_factor,inc_factor)

% Initialize Enew
Eold = 0;

% Calculate initial forces and energy
Enew = potential_energy_calc(init_coords,hp); % Initialize Enew
coords = init_coords;

% start while loop, while energy difference is larger than sigma
while abs(Enew-Eold) > sigma

    % Calculate energies and forces from current config
    Eold = potential_energy_calc(coords,hp); % Initialize Enew
    forces = force_calc(coords,hp);

    % update positions
    new_coords = coords + stepsize * forces;

    % Recalculate energy for new configuration
    Enew = potential_energy_calc(new_coords,hp);

    if (Enew < Eold)
        coords = new_coords;
        stepsize = stepsize*inc_factor;
    else
        stepsize = stepsize * dec_factor;
    end

end

min_coords = coords;
minE = potential_energy_calc(coords,hp);

end