function[KE, temp] = velocities_to_energy(velocities, mass)
    [num_particles, ~] = size(velocities);
    vel_sqsum = sum(sum(velocities.^2));
    temp = vel_sqsum/(3*num_particles - 3);
    KE = mass*vel_sqsum/2;
end