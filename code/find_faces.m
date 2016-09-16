function [dof_f_ext,dof_f] = find_faces(cell_dofs,D,Grid)
% author: Marc Hesse
% date: 7 Sep 2016
% Description:
% Find the dofs of the faces associated with a given vector of cell dof's
% and separates them into interior and exterior faces.
% The information which faces belong to which cell is in the divergence 
% matrix.
% This connectivity information is already in the discrete divergence
% matrix! Simply take the rows corresponding to the dof's of the cells
% and then sum them up. The dof's of the faces are all the non-zero
% columns. Summing all the rows eliminated the interior faces and the left
% over non-zero columns are the exerior faces.
% 
% Input:
% cell_dofs = column vector of cells in the active subdomain
% D         = discrete divergence matrix
% Grid      = structure containing pertinent information about the grid
%
% Output:
% dof_f     = column vector containing the dof's of all face in and on 
%             the exterior of the active subdomain
% dof_f_ext = 
% Select appropriate rows of the divergence matrix
DD = D(cell_dofs,:);

%% Exterior faces
% Sum columns of DD. In the columns corresponding to internal faces the
% entries will cancel. Therefore, only the sum of the columns corresponding
% to exterior faces will be non-zero.
dof_f_ext = Grid.dof_f(abs(sum(DD,1))>eps);

%% Interior faces
% Sum absolute values of the columns of DD. This way the colums
% corresponding to interior faces will not cancel so that the sum of the columns
% corresonding of all faces associated with cell dofs will be non-zero. 
% The interior cells are identified by comparison with the dofs of the
% exterior cells.
% dof_f_int = Grid.dof_f(sum(abs(DD),1)>eps);
% dof_f_int = setdiff(dof_f_int,dof_f_ext);

%% All subdomain faces
% Turns out we don't need the interior faces separately just need all the
% faces in the subdomain to restrict the operators and the exterior ones
% for the boundary conditions
dof_f = Grid.dof_f(sum(abs(DD),1)>eps);