% Calculates random force for langevin simulation. multiply this by randn
% to get the random force
function Rx = random_force(k, T, m, zeta, dt, num_particles) 
R_debug = 1; % Used for debugging
sigma = sqrt(2*k*T*m*zeta/dt);
Rx = R_debug * sigma.*randn(num_particles,3);

end
