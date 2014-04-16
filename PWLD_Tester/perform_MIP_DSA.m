function out = perform_MIP_DSA(ndat, solvdat, mesh, DoF, phi, phi0)
global glob
% Put Solution in cell format if it is not already
% ------------------------------------------------
dphi = cell(ndat.numberEnergyGroups,1);
for g=1:ndat.numberEnergyGroups
    dphi{g} = phi{g} - phi0{g};
end
ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
% Loop through energy groups
for g=1:ndat.numberEnergyGroups
    mat = sparse(ndg,ndg);
    rhs = zeros(ndg,1);
    % Loop through cells
    for tcell=1:mesh.TotalCells
        matID = mesh.get_cell_material_id(tcell);
        dofs = DoF.getLocalNodes(tcell);
        cfaces = mesh.get_cell_faces(tcell);
        [M,K,~] = PWLD_volume(mesh.get_cell_verts(tcell), mesh.get_cell_faces(tcell));
        mat(dofs,dofs) = mat(dofs,dofs) + ndat.Diffusion.DiffXS(matID,g)*sparse(K) + ndat.Diffusion.AbsorbXS(matID,g)*sparse(M);
        rhs(dofs) = rhs(dofs) + ndat.Diffusion.ScatteringXS(matID,1,g,g)*M*dphi{g}(dofs);
        % Loop through faces/edges in a cell
        for ff=1:length(cfaces)
            f = cfaces(ff);
            farea = mesh.get_face_areas(f);
            fcells = mesh.get_face_cells(f);
            fflag = mesh.get_face_flag(f);
            M = get_face_mass_matrix(mesh.Dimension,mesh.get_face_area(f));
            if fflag == 0 % interior face
                % Get info for each face cell
                % There are only 2 cells (if an interior face) no matter 
                % the problem dimension.
                matids = mesh.get_cell_material_id(fcells);
                for c=1:2
                    fcnodes{c} = DoF.get_face_cell_nodes(f,c);
                    cverts{c} = mesh.get_cell_verts(fcells(c));
                    h(c) = get_orthogonal_length(mesh.Dimension,cverts{c},farea,mesh.get_cell_faces(fcell(c)));
                    D(c) = ndat.Diffusion.DiffXS(matids(c),g);
                end
                kp = get_penalty_coefficient(DoF.Degree,D,h,fflag);
            else % boundary face
                for c=1:2
                    if fcells(c) > 0
                        fcnodes = DoF.get_face_cell_nodes(f,c); nfcnodes = length(fcnodes);
                        matids = mesh.get_cell_material_id(fcells(c));
                        cverts = mesh.get_cell_verts(fcells(c));
                        h = get_orthogonal_length(mesh.Dimension,cverts,farea,mesh.get_cell_faces(fcell(c)));
                        D = ndat.Diffusion.DiffXS(matids,g);
                    end
                end
                kp = get_penalty_coefficient(DoF.Degree,D,h,fflag);
                % Switch between boundary types
                switch(ndat.Diffusion.BCFlags(g,fflag))
                    case(glob.Dirichlet)
                        rhs(fcnodes) = rhs(fcnodes) + kp*ndat.Diffusion.BCVals(g,fflag)*M*ones(fcnodes,1);
                    case(glob.Neumann)
                        rhs(fcnodes) = rhs(fcnodes) - ndat.Diffusion.BCVals(g,fflag)*M*ones(fcnodes,1);
                    case(glob.Robin)
                        rhs(fcnodes) = rhs(fcnodes) + 2*ndat.Diffusion.BCVals(g,fflag)*M*ones(fcnodes,1);
                end
            end
        end
    end
    out{g} = mat\rhs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_orthogonal_length(dim, verts, area, faces)
if dim == 1
    out = abs(verts(2) - verts(1));
elseif dim == 2
    nv = size(verts,1);
    L = area;
    [g,~,~] = polygeom(verts(:,1),verts(:,2));
    if nv == 3
        out = 2*g(1)/L;
    elseif nv == 4
        out = g(1)/L;
    elseif nv > 4 && mod(nv,2) == 0
        out = 4*g(1)/g(4);
    elseif nv > 4 && mod(nv,2) ~= 0
        out = 2*g(1)/g(4) + sqrt((2*g(1))/(nv*sin(2*pi/nv)));
    end
elseif dim == 3
    % not implemented yet...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_penalty_coefficient(p,D,h,eflag)
if eflag == 0
    c = 2*(1+p)^2;
    out = max(0.25, c/2*(D(1)/h(1) + D(2)/h(2)));
else
    c = 4*(1+p)^2;
    out = max(0.25, 2*c*D/h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_mass_matrix(dim,area,verts)
if dim == 2
    out = area/6*[2,1;1,2];
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
