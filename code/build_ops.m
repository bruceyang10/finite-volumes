function [D,G,I]=build_ops(Grid)
% author: Marc Hesse
% date: 09/08/2014
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% I = identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

if (Nx>1) && (Ny>1) % 2D case
    %% One dimensinal divergence
    Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
    
    %% Two-dimensional divergence
    Dy = spblkdiag(Dy,Nx); % y-component
    
    e =   ones(Ny*(Nx+1),1);
    Dx = spdiags([-e e]/Grid.dx,[0 Ny],N,(Nx+1)*Ny); % x-component
    D = [Dx Dy];
    
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax;... % boundary faces
                 Grid.dof_f_ymin; Grid.dof_f_ymax];
elseif (Nx>1) && (Ny==1)
    D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1)
    D = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax];   % boundary faces
end


%% Gradient
G = -D';
G(dof_f_bnd,:) = 0;

%% Identity
I = speye(Grid.N);

% % Boundary faces (update to use Gri.dof_f!)
% if (Nx>1) && (Ny>1) % 2D case
%     bnd_x = [1:Ny,Grid.Nfx-Ny+1:Grid.Nfx];
%     bnd_y = [1:Ny+1:Grid.Nfy-Ny+1]; bnd_y = [bnd_y bnd_y+Ny];
%     bnd = [bnd_x Grid.Nfx+bnd_y];
% elseif (Nx>1) && (Ny==1)
%     bnd = [1:Ny,Grid.Nfx-Ny+1:Grid.Nfx]; % Too complicated?
% elseif (Nx==1) && (Ny>1)
%     bnd = [1 Ny+1];
% end


    


%% Adjust divergence for different coordinate systems
if strcmp(Grid.geom,'polar1D')
%     Rf = comp_mean(Grid.xc,-1, Grid);
    Rf = spdiags(Grid.xf,0,Nx+1,Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Nx,Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'spherical1D')
    Rf = spdiags(Grid.xf.^2,0,Grid.Nx+1,Nx+1);
    Rcinv = spdiags(1./(Grid.xc.^2),0,Nx,Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'cylindrical_rz')
    % assumes: y-dir is radial direction 
    %          simplifies the assembly because grid is ordered y-first
    % The change in geometry goes into 1D matrix before Dy is reassembled
    Rf = spdiags(Grid.yf,0,Ny+1,Ny+1);
    Rcinv = spdiags(1./Grid.yc,0,Ny,Ny);
    Dy = Rcinv*spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1)*Rf;
    Dy = spblkdiag(Dy,Nx);
    D = [Dx Dy];
end
