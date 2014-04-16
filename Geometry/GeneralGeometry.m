classdef GeneralGeometry < handle
    properties (Access = public)
        Dimension
        TotalCells
        TotalFaces
        TotalInteriorFaces
        TotalBoundaryFaces
    end
    properties (Access = public)
        MeshType
        Vertices
        
        MatID
        CellVerts
        CellCenter
        CellVolume
        CellFaces
        
        Faces
        InteriorFaces
        BoundaryFaces
        FaceID
        FaceNormal
        FaceCenter
        FaceArea
        FaceCells
    end
    properties (Access = public) % AMR Variables
        RefinementLevel
        RefinementBool
        RefinementArray
    end
    methods % Constructors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GeneralGeometry (varargin)
            n = nargin;
            if n == 0
                % empty constructor -> do nothing
            else
                ttime = tic;
                obj.Dimension = varargin{1};
                % 1D Construction
                if obj.Dimension == 1
                    if ~isa(varargin{2}, 'double')
                        error('1D input needs to consist of a vector of vertex coordinates.')
                    else
                        obj = Constructor_1D(varargin{2});
                    end
                elseif obj.Dimension == 2
                    if strcmp(varargin{2}, 'Triangle') || strcmp(varargin{2}, 'triangle')...
                            || strcmp(varargin{2}, 'tri') || strcmp(varargin{2}, 'Tri')
                        [verts, cells, faces] = readTriangleMesh(varargin{3});
                        [tri_out] = process_triangle_file(cells,verts,faces);
                        obj.MeshType = 'Triangle';
                        obj.CellVerts = tri_out.CellVerts;
                        obj.TotalCells = tri_out.TotalCells;
                        obj.TotalFaces = tri_out.TotalFaces;
                        obj.Vertices = tri_out.Vertices;
                        obj.Faces = tri_out.Faces;
                        obj.MatID = tri_out.MatID;
                        obj.FaceID = tri_out.FaceID;
                        obj.CellCenter = tri_out.CellCenter;
                        obj.CellVolume = tri_out.CellVolume;
                        obj.CellFaces = tri_out.CellFaces;
                        obj.FaceCenter = tri_out.FaceCenter;
                        obj.FaceNormal = tri_out.FaceNormal;
                        obj.FaceArea = tri_out.FaceArea;
                        obj.FaceCells = tri_out.FaceCells;
                    else

                    end
                elseif obj.Dimension == 3
                    
                end
                
                % Cleanup Face Terms
                n = 1; nn=1;
                for f=1:length(obj.Faces)
                    if obj.FaceID(f) == 0
                        IF(nn) = f;
                        nn = nn + 1;
                    else
                        BF(n) = f;
                        n = n + 1;
                    end
                end
                obj.InteriorFaces = uint32(IF');
                obj.BoundaryFaces = uint32(BF');
                obj.TotalInteriorFaces = length(obj.InteriorFaces);
                obj.TotalBoundaryFaces = length(obj.BoundaryFaces);
                
                % Populate Initial ARM Variables
                obj.RefinementBool = 0;
                obj.RefinementLevel = 0;
                obj.RefinementArray = zeros(obj.TotalCells,2,'uint32');
                
                disp(['Total Geometry Generation Time:  ',num2str(toc(ttime))])
                disp(' ')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = Constructor_1D(verts)
            obj.Vertices = verts;
            obj.TotalCells = length(obj.Vertices) - 1;
            obj.CellVerts = cell(obj.TotalCells,1);
            obj.Faces = cell(obj.TotalCells+1,1);
            for c=1:obj.TotalCells
                obj.CellNodes{c} = [c,c+1];
                obj.CellVerts{c} = obj.Vertices([c,c+1]);
                obj.CellFaces{c} = [c,c+1];
                obj.Faces{c} = c;
            end
            obj.Faces{end} = obj.TotalCells + 1;
            obj.FaceArea = ones(length(obj.Faces),1);
        end
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Accessor Routines - analogous to 'get' commands in C++
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_vertices(obj)
            out = obj.Vertices;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_faces(obj)
            out = obj.Faces;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_cell_verts(obj)
            out = obj.CellVerts;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_cell_faces(obj)
            out = obj.CellFaces;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_material_id(obj, cellID)
            out = obj.MatID(cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_verts(obj, cellID)
            out = obj.Vertices(obj.CellVerts{cellID},:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_vert_indices(obj, cellID)
            out = obj.CellVerts{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_faces(obj, cellID)
            out = obj.CellFaces{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_centers(obj, cellID)
            out = obj.CellCenters(cellID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_volumes(obj, cellID)
            out = obj.CellVolumes(cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_interior_face(obj, x)
            out = obj.InteriorFaces(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_boundary_face(obj, x)
            out = obj.BoundaryFaces(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_verts(obj, faceID)
            out = obj.Vertices(obj.Faces{faceID},:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_vert_indices(obj, faceID)
            out = obj.Faces{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_cells(obj, faceID)
            out = obj.FaceCells(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_normals(obj, faceID)
            out = obj.FaceNormal(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_areas(obj, faceID)
            out = obj.FaceArea(faceID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_flags(obj, faceID)
            out = obj.FaceID(faceID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_centers(obj, faceID)
            out = obj.FaceCenter(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_mesh_type(obj)
            out = obj.MeshType;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Generation Routines (i.e. set flags)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_cell_matIDs_inside_domain(obj, val, verts, faces)
            % Function to loop through cells and determine if cell
            % center lies within shape defined by 'verts' and 'faces'.
            % 'faces' is only needed in 3 dimensions.
            % -------------------------------------------------------------
            if obj.Dimension == 1
                for c=1:obj.TotalCells
                    if obj.CellCenter(c) > verts(1) && obj.CellCenter(c) < verts(2)
                        obj.MatID(c) = val;
                    end
                end
            elseif obj.Dimension == 2
                [in] = inpoly(obj.CellCenter,verts);
                for c=1:obj.TotalCells
                    if logical(in(c))
                        obj.MatID(c) = val;
                    end
                end
            elseif obj.Dimension == 3
                FV.verts = verts; FV.faces = faces;
                bools = inpolyhedron(FV,obj.CellCenter);
                for c=1:obj.TotalCells
                    if logical(bools(c))
                        obj.MatID(c) = val;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_face_flag_on_surface(obj, val, verts)
            if obj.Dimension == 1
                
            elseif obj.Dimension == 2
                for f=1:obj.TotalBoundaryFaces
                    
                end
            elseif obj.Dimension == 3
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Unstructured Routines - this is not tested for AMR...
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clear_refinement_array(obj)
            obj.RefinementArray = zeros(obj.TotalCells,2,'uint32');
            obj.RefinementBool = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_refinement_flag(obj, cellID, flag, val)
            obj.RefinementBool = 1;
            if nargin < 3
                obj.RefinementArray(cellID,:) = [flag, 1];
            else
                obj.RefinementArray(cellID,:) = [flag, val];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_mesh(obj)
            for c=1:obj.TotalCells
                if obj.RefinementArray(c,1) ~= 0
                    refine_cell(c);
                end
            end
            obj.RefinementLevel = obj.RefinementLevel + 1;
            obj.clear_refinement_array();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_cell(obj, cellID)
            if obj.Dimension == 1
                obj.refine_cell_1D(obj, cellID);
            elseif obj.Dimension == 2
                obj.refine_cell_2D(obj, cellID);
            elseif obj.Dimension == 3
                obj.refine_cell_3D(obj, cellID);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_cell_1D(obj, cellID)
            val = obj.RefinementArray(cellID,2);
            
            cverts = obj.CellVerts{cellID};
            verts = obj.Vertices(cverts,:);
            dx = abs(verts(1) - verts(2));
            nc = 2^val;
            ddx = dx/nc;
            newc = nc-1;
            for r=1:nc
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_cell_2D(obj, cellID)
            flag = obj.RefinementArray(cellID,1);
            val = obj.RefinementArray(cellID,2);
            switch(flag)
                case(1) % Add new point in cell center
                    verts = obj.Vertices(obj.CellVerts{cellID},:);
                    rcenter = mean(verts);
                    faces = obj.CellFaces{cellID};
                    newcells = length(faces)-1;
                    newfaces = newcells;
                    % Loop through nverts-1
                    for i=1:newcells
                        nverts = [verts([i,i+1],:);rcenter];
                        
                    end
                    % Last cell addition
                    
                    obj.TotalCells = obj.TotalCells + newcells;
                case(2)
                    verts = obj.Vertices(obj.CellVerts{cellID},:);
                    if size(verts,1) == 3
                        ncells = 4;
                    elseif size(verts,1) == 4
                        ncells = 4;
                    else
                        error('This refinement option only works for triangles and quads.')
                    end
                otherwise
                    
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_cell_3D(obj, cellID)
            flag = obj.RefinementArray(cellID,1);
            val = obj.RefinementArray(cellID,2);
            switch(flag)
                case(1)
                    
                case(2)
                    
                otherwise
                    
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Accessory Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = get_face_info(obj, verts)
    
end


