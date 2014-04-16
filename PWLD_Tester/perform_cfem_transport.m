function out = perform_cfem_transport(ndat, solvdat, mesh, DoF, x)

% Loop through Quadrature Directions
% ----------------------------------
for m=1:ndat.Transport.NumberQuadratureDirections
    [L, rhs] = get_global_matrices(angDir, ndat, mesh, DoF, x);
    y = L\rhs;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, rhs] = get_global_matrices(angNum, ndat, mesh, DoF, x)
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
ndim = mesh.Dimension;
angDir = ndat.Transport.AngularDirections(angNum);
angWgt = ndat.Transport.AngularWeights(angNum);
% Allocate Memory
L = sparse(ndof,ndof);
rhs = zeros(ndof,1);
% Loop through cells
for cell=1:mesh.TotalCells
    matID = mesh.get_cell_material_id(cell);
    [M,~,G] = PWLD_volume(mesh.get_cell_verts(cell), mesh.get_cell_faces(cell));
    for g=1:ng
        sg = DoF.getShiftedLocalNodes((g-1)*ndg,cell);
        L(sg,sg) = sparse(L(sg,sg)) + ndat.Transport.TotalXS(matID,g)*M;
        for d=1:ndim
            L(sg,sg) = sparse(L(sg,sg)) - angDir(d)*G{d};
        end
        rhs(sg) = rhs(sg) + ndat.Transport.ExtSource(matID,g)*M*ones(length(sg),1);
        for gg=1:ng
            sgg = DoF.getShiftedLocalNodes((gg-1)*ndg,cell);
            rhs(sg) = rhs(sg) + ndat.FissSpec(matID,g)/ndat.keff*ndat.FissionXS(matID,gg)/(4*pi)*M*x(sgg);
            for m=1:ndat.Transport.fluxMoments;
                mm = m-1;
                sxs = (2*mm+1)/(4*pi)*ndat.Transport.ScatteringXS(matID,gg,g,m);
                rhs(sg) = rhs(sg) + sxs*M*x(sgg);
%                 for mmm=-m:m
%                     mmmm = mmm + (2*m+1);
%                 end
            end
        end
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.get_boundary_face(f);
    fnorm = mesh.get_face_normals(face);
    if dot(angDir,fnorm) < 0
        flag = mesh.get_face_flags(face);
        tflag = ndat.Transport.BCFlags(flag);
        switch(tflag)
            case(out.Vacuum)
                
            case(out.Reclecting)
                
            case(out.IncidentIsotropic)
                
            case(out.IncidentCurrent)
                
        end
        for g=1:ng
            fcnodes = DoF.getShiftedFaceCellNodes((g-1)*ndg,face,1);
            fflag = ndat.Transport.BCFlags(flag);
        end
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_source(m)

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%