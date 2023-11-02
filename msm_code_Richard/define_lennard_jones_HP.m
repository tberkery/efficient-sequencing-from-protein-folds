function[sigma, epsilon]=define_lennard_jones_HP(mol1, mol2)
    sigma = 1;
    epsilon = 2/3;
    if (mol1 == 'H' && mol2 == 'H')
        epsilon = 1;
    end
end