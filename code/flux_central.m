function [A] = flux_central(q,Grid)
% q =  volumetric flux on faces
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces

if Ny == 1 && Nz == 1 % 1D
    %% One dimensinal

    A = spdiags([q(1:Nx),q(2:Nx+1)]/2,[-1 0],Grid.Nx+1,Grid.Nx);
else
    error('2D central flux is not implementd.')
end