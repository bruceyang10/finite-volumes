function [u] = solve_lbvp(L,f,B,g,N)
% author: Marc Hesse
% date: 26 Sept 2014
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in \Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
%
% Input:
% L = matrix representing the discretized linear operator of size N by N, 
%     where N is the number of degrees of fredom
% f = column vector representing the discretized r.h.s. and contributions
%     due non-homogeneous Neumann BC's of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC's of
%     size Nc by N
% g = column vector representing the non-homogeneous Dirichlet BC's of size
%     Nc by 1.
% N = matrix representing a orthonormal basis for the null-space of B and
%     of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> dof_dir = 1;
% >> B = I(dof_dir,:); g = 1; 
% >> N = I; N(:,dof_dir) = [];
% >> h = solve_lbvp(L,fs,B,g,N);

if isempty(B) % no constraints
    u = L\f;
else
    %u0 = zeros(length(f)); 
    up = spalloc(length(f),1,length(g));
    up = B'*(B*B'\g);
    u0 = N*(N'*L*N\(N'*(f-L*up)));
    u = u0 + up;
end