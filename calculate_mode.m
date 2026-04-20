function psi = calculate_mode(Nx,Ny,Nz,Lx,Ly,Lz,r)
    % Calculate a single mode
    % Assumes perfect pressure-release boundaries at all walls and surface
    
    % Inputs:
    % Nx, Ny, Nz: Mode number in each dimension
    % Lx, Ly, Lz: Tank dimensions (m)
    % r = [rx; ry; rz]: vector position (m) [3x1] or [1x3]

    % Outputs:
    % psi: tank eigenvector
    
    psi = sin(Nx.*pi.*r(1)./Lx).*sin(Ny.*pi.*r(2)./Ly).*sin(Nz.*pi.*r(3)./Lz);
end