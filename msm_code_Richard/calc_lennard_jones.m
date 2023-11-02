function[lennard_jones]=calc_lennard_jones(distances, sigma, epsilon)
    lennard_jones = 0;
    for distance = distances
        lennard_jones = lennard_jones + 4 * epsilon * ...
            ((sigma / distance)^12 - (sigma / distance)^6);
    end
    
end