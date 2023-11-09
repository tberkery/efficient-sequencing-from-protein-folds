% Calculate forces for MD, returns a three item row vector containing
% forces for i and j
% Requires the function files MIC and lennard_force
function forces = force_calc(coords, is_hydrophil)
% Initialize matrix of the same size as the positions
[num_particles, dim] = size(coords);
forces = zeros(num_particles,dim);
% Initalize variables for lennard Jones
sigma = 1;

% Array that shows whether the particle is hydrophilic or hydrophobic
%is_hydrophil = [true;false;false;true;true];

%% Calculate total config force due to 1)lennard jones 2)harmonic bonded 
for i = 1:num_particles
  for j = (i+1):num_particles

      % Calculate forces. Fij is a 3D vector of the force exerted by j on i
      if (is_hydrophil(i) && is_hydrophil(j))
          epsilon = 1;
      else
          epsilon = 2/3;
      end

      % Calculate bonded force for adjacent particles and LJ force for
      % non-adjacent
      if (i+1) == j
          Fij = harmonic_bonded_forces(coords(i,:), coords(j,:), 20, 1);
      else
          Fij = lennard_force_noMIC(coords(i,:), coords(j,:), epsilon, sigma);
      end
      
      forces(i,:) = forces(i,:) + Fij;
      forces(j,:) = forces(j,:) - Fij;

  end
end


end