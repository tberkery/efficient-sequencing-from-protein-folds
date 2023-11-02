function[coords_minE, F_minE, E_minE]=steepest_descent_HP(hpChain, chain_coords, step_size, step_factor, tolerance)
    % curr coords, force, energy
    [npart, ndimensions] = size(chain_coords);
    numid = zeros(numel(hpChain), 1);
    for i = 1:numel(numid)
        if (hpChain{i} == 'P')
            numid(i) = 1;
        end
    end
    kh = 20;
    le = 1;
    sig = 1;
    coords_curr = chain_coords;
    [F_curr, E_curr] = calc_en_f(coords_curr, numid, npart, sig, kh, le);
    E_minE = E_curr;
    coords_minE = coords_curr;
    F_minE = F_curr;
    
    coords_new = coords_curr + step_size*F_curr;
    [F_new, E_new] = calc_en_f(coords_new, numid, npart, sig, kh, le);
    while(abs(E_new - E_curr) > tolerance)
        % minima along force (the negative gradient of Energy)
        if (E_new > E_curr)
            step_size = step_size*step_factor;
        elseif (E_new < E_curr)
            step_size = step_size/step_factor;
            
            % update
            coords_curr = coords_new;
            F_curr = F_new;
            E_curr = E_new;
            
            % update min
            if (E_new < E_minE)
                E_minE = E_new;
                coords_minE = coords_new;
                F_minE = F_new;
            end
        end
        % new coords, force, energy
        coords_new = coords_curr + step_size*F_curr;
        [F_new, E_new] = calc_en_f(coords_new, numid, npart, sig, kh, le);
    end
end