function refine_mesh(mesh, DoF, pdata, x)
if ~iscell(x)
    n = size(x,2); y=cell(n,1);
    for i=1:n
        y{i} = x(:,i);
    end
else
    y = x;
end
clear x
% Determine Mesh Refinement Based on FEM Type
%   CFEM = Gradient Jumps on Edges
%   DFEM = Solution Jumps on Edges
% -------------------------------------------
if DoF.FEMType == 1
    % Loop through faces in mesh
    % --------------------------
    for IF=1:mesh.TotalInteriorFaces
        face = mesh.get_interior_face(IF);
        cells = mesh.get_face_cells(face);
        fnodes = DoF.getFaceDoFs(face);
        for c=1:length(y)
            
        end
    end
elseif DoF.FEMType == 2
    err = zeros(mesh.TotalCells,1);
    % Loop through cells in mesh
    % Estimate error for each cell
    % ----------------------------
    for c=1:mesh.TotalCells
        faces = mesh.get_cell_faces(c);
        err(c) = DFEM_error_measure();
    end
    maxerr = max(err);
    % Loop through cells in mesh
    % Apply refinement
    % --------------------------
    for c=1:mesh.TotalCells
        if err(c) >= pdata.refinementTolerance*maxerr
            mesh.set_refinement_flag(c,pdata.refinementType,pdata.refinementSplits);
        end
    end
end

mesh.refine_mesh();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = CFEM_error_measure()

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = DFEM_error_measure()

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%