function [dof_bnd_in,dof_bnd_out] = find_bnd_cells(dof_in,dof_out,dof_f_bnd,D,Grid)

DD = D(:,dof_f_bnd);
dof_bnd_in = Grid.dof(abs(sum(DD,2))>eps); % these are all cells. Using dof_bnd_in for tmp storage

if length(dof_in) > length(dof_out)
    % Cheaper to compare with exterior cells
    dof_bnd_in = Grid.dof(abs(sum(DD,2))>eps); 
    % Note, these are all cells! Using dof_bnd_in for temporary storage
    [dof_bnd_out,i_bnd_out] = intersect(dof_bnd_in,dof_out); % find exterior cells
    dof_bnd_in(i_bnd_out) = []; % delete exterior cells
else
    % Cheaper to compare with interior cells
    dof_bnd_out = Grid.dof(abs(sum(DD,2))>eps); 
    % Note, these are all cells! Using dof_bnd_out for temporary storage
    [dof_bnd_in,i_bnd_in] = intersect(dof_bnd_out,dof_in); % find interior cells
    dof_bnd_out(i_bnd_in) = []; % delete interior cells
end


