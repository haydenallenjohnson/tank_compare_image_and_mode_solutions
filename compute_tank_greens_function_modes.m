function g_omega = compute_tank_greens_function_modes(r_source,r_receiver,omega,Lx,Ly,Lz,c,N_max,damping)
    
    % Calculate the Green's function for the tank using a modal sum.

    % Inputs:
    % r_source: Vector position of the source (m) [3x1]
    % r_receiver: Vector position of the receiver (m) [3x1]
    % omega: Angular frequencies at which to compute G (radians) [Nx1]
    % Lx, Ly, Lz: Dimensions of the tank in each coordinate (m)
    % c: Sound speed (m/s)
    % N_max: Highest mode number to include (same for all dimensions)
    % damping: Constant damping applied which simulates gradual loss of
    % energy, evenly across frequencies, proportional to distance
    % travelled
   
    % Outputs:
    % g_tank: [Nx1] Green's function for propation in the tank between the 
    % source and receiver positions, as a function of the angular 
    % frequencies specified by omega.

    % Hayden Johnson 2026-04-13
        
    % compute mode normalization factor
    Lambda = Lx*Ly*Lz./8;
    L = [Lx Ly Lz];

    % create output array
    g_omega = zeros(size(omega));
    n = length(omega);
    
    % iterate over the first n/2+1 frequencies
    parfor i = 1:n/2+1
    
        k_omega = omega(i)./c + 1i*damping;

        % add modes together
        s = 0;
        for Nx = 1:N_max
            for Ny = 1:N_max
                for Nz = 1:N_max
                    N = [Nx Ny Nz]
                    k_mode_squared = sum((pi.*N./L).^2);
                    s = s + calculate_mode(Nx,Ny,Nz,Lx,Ly,Lz,r_source).*calculate_mode(Nx,Ny,Nz,Lx,Ly,Lz,r_receiver)./(Lambda.*(k_mode_squared-k_omega.^2));
                end
            end
        end
        
        % write output value
        g_omega(i) = s;
    end
    
    % the rest of the frequencies are complex conjugates of the first half
    g_omega(n/2+2:end) = flip(g_omega(2:n/2))';
end