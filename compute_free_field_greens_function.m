function g_free = compute_free_field_greens_function(r1,r2,omega,c,damping)
    
    % Compute the free-field green's function, in the frequency domain, for
    % a point source and receiver and positions r1 and r2. (Note that r1
    % and r2 are interchangeable)
    
    % Inputs:
    % r1: Vector position of the source (m) [3x1]
    % r2: Vector position of the receiver (m) [3x1]
    % omega: Angular frequencies at which to compute g (radians) [Nx1]
    % c: Sound speed (m/s)
    % damping: Constant damping applied which simulates gradual loss of
    % energy, evenly across frequencies, proportional to distance
    % travelled
    
    % Outputs:
    % g_free: [Nx1] Green's function for propation between the source and 
    % receiver positions, as a function of the angular frequencies 
    % specified by omega.
    
    % Written by Hayden Johnson, 2024-03-11
    
    %----------------------------------------------------------------------
    
    % initialize array
    N = length(omega);
    g_free = zeros(N,1);
    
    % compute separation distance
    r = sqrt(sum((r2-r1).^2));

    % compute wavenumber
    k = omega(1:N/2+1)./c + 1i*damping;
    
    % compute first N/2+1 values of the greens function
    g_free(1:N/2+1) = exp(-1i*r*k)./(4*pi*r);
    % added factor of 1/(4*pi) 2026-04-13 to agree with modes
    
    % specify the rest of the values as the complex conjugates of the first
    % half
    g_free(N/2+2:end) = conj(flip(g_free(2:N/2))); 
end