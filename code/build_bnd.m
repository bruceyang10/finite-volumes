function [B,N,fn] = build_bnd(Param,Grid)
% author: Marc Hesse
% date: 06/09/2015
% Description:
% This function computes the operators and r.h.s vectors for both Dirichlet
% and Neumann boundary conditions.
%
% Input:
% Grid = structure containing all pertinent information about the grid.
% Param = structure containing all information about the physical problem
%         in particular this function needs the fields
%         Param.dof_dir = Nc by 1 column vector containing 
%                         the dof's of the Dirichlet boundary.
%         Param.dof_neu = N by 1 column vector containing 
%                         the dof's of the Neumann boundary.
%         Param.qb      = column vector of prescribed fluxes on Neuman bnd.
%
% Output:
% B = Nc by N matrix of the Dirichlet constraints
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B
% fn = N by 1 r.h.s. vector of Neuman contributions
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir = Grid.dof_xmin;  % identify Dirichlet bnd
% >> Param.dof_neu = Grid.dof_xmax;  % identify Neumann bnd
% >> Param.qb = 1;                   % set bnd flux
% >> [B,N,fn] = build_bnd(Param,Grid);

%% Check input format
if isrow(Param.dof_dir)   && length(Param.dof_dir)>1;   error('Param.dof_dir is not a column vector'); end
if isrow(Param.dof_neu)   && length(Param.dof_neu)>1;   error('Param.dof_neu is not a column vector'); end
if isrow(Param.dof_f_dir) && length(Param.dof_f_dir)>1; error('Param.dof_f_dir is a not column vector'); end
if isrow(Param.dof_f_neu) && length(Param.dof_f_neu)>1; error('Param.dof_f_neu is a not column vector'); end
if isfield(Param,'qb') && isrow(Param.qb) && length(Param.qb)>1;        error('Param.qb is not a column vector'); end

%% Dirichlet boundary conditions
if isempty(Param.dof_dir)
    B = [];
    N = [];
else
    N = speye(Grid.N); % use N as temp storage for identity
    B = N(Param.dof_dir,:);
    N(:,Param.dof_dir) = [];
end

%% Neumann boundary conditions
if isempty(Param.dof_neu)
    fn = spalloc(Grid.N,1,0);                     % allocate sparse zero vector
else
    fn = spalloc(Grid.N,1,length(Param.dof_neu)); % allocate sparse vector
    fn(Param.dof_neu) = Param.qb.*Grid.A(Param.dof_f_neu)./Grid.V(Param.dof_neu);
end

