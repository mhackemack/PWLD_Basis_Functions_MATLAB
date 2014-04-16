function x = perform_cfem_diffusion(ndat, solvdat, mesh, DoF, x)

ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end

[A,rhs] = get_global_matrices(x, ndat, mesh, DoF);
x = A\rhs;
x = cleanup_small_vals(x);
x = vector_to_cell(x,DoF);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat,rhs] = get_global_matrices(v, ndat, mesh, DoF)
global glob
n = length(v);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;
% Allocate Memory
rhs = zeros(n,1);
mat = sparse(n,n);
% Loop through cells
for cell=1:mesh.TotalCells
    matID = mesh.get_cell_material_id(cell);
    [M,K,~] = PWLD_volume(mesh.get_cell_verts(cell), mesh.get_cell_faces(cell));
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = DoF.getShiftedLocalNodes((g-1)*ndg,cell);
        mat(sg,sg) = sparse(mat(sg,sg)) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M;
        rhs(sg) = rhs(sg) + ndat.Diffusion.ExtSource(matID,g)*M*ones(length(sg),1);
        for gg=1:ng
            sgg = DoF.getShiftedLocalNodes((gg-1)*ndg,cell);
            rhs(sg) = rhs(sg) + ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*M*v(sgg);
            rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,gg,g,1)*M*v(sgg);
        end
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.get_boundary_face(f);
    flag = mesh.get_face_flags(face);
    % Loop through energy groups
    for g=1:ng
        fcnodes = DoF.getShiftedFaceCellNodes((g-1)*ndg,face,1);
        gflag = ndat.Diffusion.BCFlags(g,flag);
        switch(gflag)
            case(glob.Dirichlet)
                mat(fcnodes,:) = 0;
                for i=1:length(fcnodes)
                    mat(fcnodes(i),fcnodes(i)) = 1;
                end
                rhs(fcnodes) = ndat.Diffusion.BCVals(g,gflag)*ones(length(fcnodes),1);
            case(glob.Neumann)
                rhs(fcnodes) = rhs(fcnodes) - ndat.Diffusion.BCVals(g,gflag)*M*ones(length(fcnodes),1);
            case(glob.Robin)
                rhs(fcnodes) = rhs(fcnodes) + 2*ndat.Diffusion.BCVals(g,gflag)*M*ones(length(fcnodes),1);
        end
    end
end
return