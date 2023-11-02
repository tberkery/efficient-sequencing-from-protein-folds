% Melissa Mai
% Function calc_en_f
%
% Calculates potential energy and forces in Langevin Dynamics. Energy
% includes Lennard-Jones non-bonded and Harmonic bonded potential. Forces
% include Lennard-Jones, harmonic bonded, drag, and random forces.
%
%
% Inputs
% coord: coordinate matrix for configuration
% numid: numerical ID vector for particles
% vels: velocity matrix
% npart: number of particles
% T: temperature
% zeta: friction coefficient
% dt: time step
%
% Optional inputs
% sig: sigma, default = 1
% kh: harmonic bond force constant, default = 20
% le: bond length, default = 1
% m: mass, default = 1
% k: Boltzmann constant, default = 1
%
% Outputs 
% LJE: Lennard-Jones potential
% HE: harmonic potential
% Ftot: total force matrix
% Fc: total deterministic forces
% Flj: LJ forces
% Fh: harmonic forces
% Fr: random forces
% Fd: drag forces

function [Ftot, Etot] = calc_en_f(coords, numid,...
    npart, sig, kh, le)

    % Assign default values if not given by user
    if nargin < 8
        sig = 1;
        kh = 20;
        le = 1;
    end
    
    % Initialize variables
    LJE = 0; HE = 0;
    Flj = zeros(npart, 3);
    Fh = zeros(npart, 3);
        
    % Calculate pairwise interactions
    for i = 1:npart
        parti = coords(i,:); 
        for j = i+1:npart
            partj = coords(j,:);
            
            % Distance between particles
            xij = parti(1) - partj(1);
            yij = parti(2) - partj(2);
            zij = parti(3) - partj(3);
            
            rij2 = xij^2 + yij^2 + zij^2;
            
            % Determine interaction type and corresponding epsilon value
            % Check for H-H interaction
            if numid(i) + numid(j) > 0 % is not H-H interaction
                eps = 2/3;
            else % is H-H interaction
                eps = 1;
            end            
            
            % Energy Calculation
            if abs(j-i) > 1 % if particles are not adjacent to each other
                r2i = (sig^2)/rij2;
                r6i = r2i^3;
                
                % LJ Potential calculation
                LJE = LJE + 4*eps*(r6i^2 - r6i);
                
                % LJ force calculation
                flj = 48 * eps * r2i * (r6i^2 - 0.5*r6i);
                Flj(i,:) = Flj(i,:) + flj*[xij yij zij];
                Flj(j,:) = Flj(j,:) - flj*[xij yij zij];
            else % if particles are adjacent to each other
                rij = sqrt(rij2);
                
                % Harmonic potential calculation
                HE = HE + 0.5*kh*(rij - le)^2;
                
                % Harmonic force calculation
                fh = -kh*(rij-le);
                Fh(i,:) = Fh(i,:) + fh*[xij yij zij]./rij;
                Fh(j,:) = Fh(j,:) - fh*[xij yij zij]./rij;
            end
        end
    end
    
    
    % Sum deterministic forces
    Fc = Flj + Fh;
    
    % Total forces
    Ftot = Fc;
    Etot = LJE + HE;
    