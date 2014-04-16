function x = perform_dfem_diffusion(ndat, solvdat, mesh, DoF, x)

ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end

[L,rhs] = get_global_matrices(ndat, mesh, DoF, x);
x = L\rhs;
x = cleanup_small_vals(x);
x = vector_to_cell(x,DoF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_global_matrices(ndat, mesh, DoF, x)
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
ndim = mesh.Dimension;
% Allocate Memory
L = sparse(ndof,ndof);
rhs = zeros(ndof,1);
% Loop through cells
for tcell=1:mesh.TotalCells
    ncnodes = length(DoF.getLocalNodes(tcell));
    cfaces = mesh.get_cell_faces(tcell);
    matID = mesh.get_cell_material_id(tcell);
    [M,K,~] = PWLD_volume(mesh.get_cell_verts(tcell), mesh.get_cell_faces(tcell));
    [Mf,Gf] = PWLD_surface_ind(mesh.get_cell_verts(tcell), mesh.get_cell_faces(tcell));
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = DoF.getShiftedLocalNodes((g-1)*ndg,tcell);
        L(sg,sg) = sparse(L(sg,sg) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M);
        rhs(sg) = rhs(sg) + ndat.Diffusion.ExtSource(matID,g)*M*ones(length(sg),1);
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = DoF.getShiftedLocalNodes((gg-1)*ndg,tcell);
            rhs(sg) = rhs(sg) + ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*M*x(sgg);
            rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,gg,g,1)*M*x(sgg);
        end
    end
    % Loop through faces
    for ff=1:length(cfaces)
        D = []; fcnodes = [];
        f = cfaces(f);
        fcells = mesh.get_face_cells(f);
        farea = mesh.get_face_areas(f);
        fnorm = mesh.get_face_normals(f);
        fflag = mesh.get_face_flag(f);
        tG = norm_gradient_mat(ndim, ncnodes, fnorm, Gf{ff});
        cmult = get_cell_mult(tcell,fcells);
        if flag == 0 % Interior Face
            matids = mesh.get_cell_material_id(fcells);
            for c=1:2
                fcnodes{c} = DoF.get_face_cell_nodes(f,c);
                cverts{c} = mesh.get_cell_verts(fcells(c));
                h(c) = get_orthogonal_length(mesh.Dimension,cverts{c},farea,mesh.get_cell_faces(fcell(c)));
                for g=1:ndat.numberEnergyGroups
                    D(g,c) = ndat.Diffusion.DiffXS(matids(c),g);
                end
            end
            % Apply interior terms
            for g=1:ndat.numberEnergyGroups
                kp = get_penalty_coefficient(DoF.Degree,D(g,:),h,fflag);
                sg = DoF.getShiftedFaceCellNodes((g-1)*ndg,f,1); nfcnodes = length(sg);
                
            end
        else % Exterior Face
            for c=1:2
                if fcells(c) > 0
                    fcnodes = DoF.get_face_cell_nodes(f,c); 
                    matids = mesh.get_cell_material_id(fcells(c));
                    cverts = mesh.get_cell_verts(fcells(c));
                    h = get_orthogonal_length(mesh.Dimension,cverts,farea,mesh.get_cell_faces(fcell(c)));
                    for g=1:ndat.numberEnergyGroups
                        D(g) = ndat.Diffusion.DiffXS(matids,g);
                    end
                    break
                end
                % Apply boundary terms
                for g=1:ndat.numberEnergyGroups
                    kp = get_penalty_coefficient(DoF.Degree,D(g),h,fflag);
                    sg = DoF.getShiftedFaceCellNodes((g-1)*ndg,f,1); nfcnodes = length(sg);
                    % Switch between boundary types
                    switch(ndat.Diffusion.BCFlags(g,fflag))
                        case(glob.Dirichlet)
                            rhs(sg) = rhs(sg) + kp*ndat.Diffusion.BCVals(g,fflag)*Mf{ff}*ones(nfcnodes,1);
                            rhs(sg) = rhs(sg) - ndat.Diffusion.BCVals(g,fflag)*D(g)*tG'*ones(nfcnodes,1);
                            mat(sg,sg) = sparse(mat(sg,sg) + kp*Mf{ff} - D(g)*(tG + tG'));
                        case(glob.Neumann)
                            rhs(sg) = rhs(sg) - ndat.Diffusion.BCVals(g,fflag)*Mf{ff}*ones(nfcnodes,1);
                        case(glob.Robin)
                            rhs(sg) = rhs(sg) + 2*ndat.Diffusion.BCVals(g,fflag)*Mf{ff}*ones(nfcnodes,1);
                            mat(sg,sg) = sparse(mat(sg,sg) + (1/2)*Mf{ff});
                    end
                end
            end
        end
    end
end
return
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
return
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
function out = norm_gradient_mat(dim, nv, fnorm, mat)
out = zeros(nv,nv);
for d=1:dim
    out = out + fnorm(d)*mat{d};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_cell_mult(tcell, fcells)
if tcell == fcells(1)
    out = -1;
else
    out = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_conforming_mat(mat)
out = flipud(mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%