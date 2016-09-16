function [dof_in,dof_out] = find_poly_cells(x_poly,y_poly,Grid)

[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
in = inpolygon(Xc(:),Yc(:),x_poly,y_poly);
dof_in  = Grid.dof(in==1);
dof_out = Grid.dof(in==0);

