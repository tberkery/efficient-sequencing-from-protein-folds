function [forces, PE] = calc_forces_configuration(hpChain, coords, k, le)
    % Initialize variables
    [num_particles, num_dimensions] = size(coords);
    PE = 0;
    forces = zeros(num_particles, num_dimensions);
    
    % Iterate through each pair of particles
    for i = 1:num_particles
        for j = i+1:num_particles
            % Distance between atoms in the x, y, and z directions
            r_ij = calc_distance_vector(coords(i, :), coords(j, :));
            r2 = sum(r_ij.^2);
            
            if (j == i + 1) % harmonic bond
                % force
                f = k*(le - sqrt(r2));
                forces(i,:) = forces(i,:) + f*r_ij;
                forces(j,:) = forces(j,:) - f*r_ij;
                % PE
                PE = PE + calc_potential_harmonicBonded(r2, k, le);
            else % lennard jones
                [sigma, epsilon] = define_lennard_jones_HP(hpChain{i}, hpChain{j});

                % Squared distance in 3 dimensions
                r2i = (sigma^2)/r2;
                r6i = r2i^3;

                % Force field determined from F = -dV/dr
                f = 48 * epsilon * r2i * (r6i^2 - 0.5*r6i);

                % Add the force on particles i and j from each other to
                % the other forces on them from other particles. The force
                % calculated above is multiplied by dx, dy, or dz to scale
                % the force in each direction.
                forces(i,:) = forces(i,:) + f*r_ij;
                forces(j,:) = forces(j,:) - f*r_ij;

                % Add potential energy of this pair to the total potential
                % energy of the system
                PE = PE + (r6i^2 - r6i);
            end
        end
    end
end