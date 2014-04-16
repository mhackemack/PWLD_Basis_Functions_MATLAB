function plot_2D_solution(mesh, DoF, x)

mtype = mesh.get_mesh_type();
ftype = DoF.FEMName;

if strcmp(ftype, 'CFEM')
    if strcmp(mtype, 'Triangle')
        plot_triangle_cfem_solution(mesh,x);
    else
        
    end
else
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_triangle_cfem_solution(mesh,x)
verts = mesh.get_all_vertices();
cells = mesh.get_all_cell_verts();
tcells = convert_cell_to_elements(cells);
if iscell(x), x = x{1}; end
trisurf(tcells,verts(:,1),verts(:,2),x)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = convert_cell_to_elements(cells)
m = length(cells{1});
n = length(cells);
out = zeros(n,m);
for i=1:n
    out(i,:) = cells{i};
end
return