% Calculate potential energy for MD, LD or other coarse grain models
% Combines lennard jones potential

function PE = potential_energy_calc(coords,hp)

% Initialize matrix of the same size as the positions
[num_particles, ~] = size(coords);
% Initalize variables for lennard Jones
sigma = 1;

%Iterate through each pair and calculate lennard jones potential
PE = 0;
for i = 1:num_particles
  x1 = coords(i,1);
  y1 = coords(i,2);
  z1 = coords(i,3);
  for j = (i+2):num_particles
      x2 = coords(j,1);
      y2 = coords(j,2);
      z2 = coords(j,3);
      if (hp(i)=='H' && hp(j)=='H')
          epsilon = 1;
      else
          epsilon = 2/3;
      end
      PE = PE + lennard_jones_noMIC(x1,y1,z1,x2,y2,z2,epsilon,sigma);
  end
end

% Harmonic bonded PE
% Iterate through each bonded pair and calculate bonded potential
bonded_PE = 0;
for i = 1:(num_particles-1)
    bonded_PE = bonded_PE + harmonic_bonded_PE(coords(i,:), coords(i+1,:),20,1);
end
PE = PE + bonded_PE;