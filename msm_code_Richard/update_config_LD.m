function [newcoords, newvels, newforces, newenpot, newenkin]...
    = update_config_LD(hpChain, coords, vels, forces, dt, temp, mass, friction, k, le)
    
    % Initialize variables
    [num_particles, num_dimensions] = size(coords);
    
    % Verlet position update
    newcoords = coords + dt*vels + 0.5*(dt^2)*forces/mass;
    
    % Calculate force field and potential energy for new configuration
    %%%%%[newforces, newenpot] = calc_forces_LD(hpChain, coords, vels, dt, temp, mass, friction, k, le);
    [forces_config, PE_config] = calc_forces_configuration(hpChain, newcoords, k, le);
    forces_random = randn(num_particles, num_dimensions) * sqrt(2*temp*mass*friction/dt);
    
    % Velocity verlet with new forces
    delta = dt/(2*mass);
    newvels = 1/(1+(dt*friction/2)) * ...
        (vels + delta*forces + delta*(forces_config + forces_random));
    
    forces_dissipative = friction*mass.*newvels;
    newforces = forces_config + forces_random - forces_dissipative;
    newenpot = PE_config;
    
    % Evaluate new temperature and kinetic energy
    [newenkin, ~] = velocities_to_energy(newvels, mass);
end