function [potentials_harmonicBonded] = calc_potential_harmonicBonded(...
    distances, k, le)
    potentials_harmonicBonded = 0.5*k*((distances - le).^2);
end