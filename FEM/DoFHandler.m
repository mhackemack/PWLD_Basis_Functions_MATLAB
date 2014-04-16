classdef DoFHandler < handle
    properties (Access = public)
        TotalDoFs
        Degree
        FEMType % 1 = CFEM, 2 = DFEM
        FEMName
    end
    properties (Access = public)
        Dimension
        GeometryType
        
        TotalCells
        TotalFaces
        
        FaceCellNodes
        CellFaceNodes
        ConnectivityArray
        NodeLocations
        ConstraintArray
    end
    properties (Access = private)
        
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                           Constructor
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = DoFHandler (varargin)
            n = nargin;
            if n == 0
                % empty constructor -> do nothing
            elseif n == 3
                % Read in Geometry and Finite Element Order
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                d = varargin{1};
                obj.Dimension = d.Dimension;
                obj.TotalCells = d.TotalCells;
                obj.TotalFaces = d.TotalFaces;
                obj.Degree = varargin{2};
                obj.CellFaceNodes = cell(obj.TotalCells,1);
                obj.FaceCellNodes = cell(obj.TotalFaces,2);
                obj.ConnectivityArray = cell(obj.TotalCells,1);
                if strcmp(lower(varargin{3}),'cfem')
                    obj.FEMType = 1;
                    obj.FEMName = 'CFEM';
                elseif strcmp(lower(varargin{3}),'dfem')
                    obj.FEMType = 2;
                    obj.FEMName = 'DFEM';
                else
                    error('Unsure of FEM Type.')
                end
                
                clear varargin
                                
                % Get Geometry Information and Build DoF Arrays
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.Dimension == 1
                    obj.GeometryType = 'Quads';
                elseif obj.Dimension == 2
                    if strcmp(d.get_mesh_type(),'Triangle')
                        obj.GeometryType = 'Triangles';
                    elseif strcmp(d.get_mesh_type(),'Quadrilateral')
                        obj.GeometryType = 'Quads';
                    else
                        obj.GeometryType = 'Polygons';
                    end
                elseif obj.Dimension == 3
                    if strcmp(d.get_mesh_type(),'Tetrahedra')
                        obj.GeometryType = 'Triangles';
                    elseif strcmp(d.get_mesh_type(),'Hexahedra')
                        obj.GeometryType = 'Quads';
                    else
                        obj.GeometryType = 'Polygons';
                    end
                else
                    error('You specified more than 3 spatial dimensions...')
                end
                
                if obj.Dimension == 1
                    obj = generate1DDoFs(obj, d);
                elseif obj.Dimension == 2
                    obj = generate2DDoFs(obj, d);
                elseif obj.Dimension == 3
                    obj = generate3DDoFs(obj, d);
                end
            else
                error('Wrong number of constructor arguments...')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                         Manipulators
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                            Accessors
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getTotalCells (obj)
            out = obj.TotalCells;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getDimension (obj)
            out = obj.Dimension;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getDegree (obj)
            out = obj.Degree;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getGeometryType (obj)
            out = obj.GeometryType;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getConstraintArray (obj)
            out = obj.ConstraintArray;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getNodeLocations (obj)
            out = obj.NodeLocations;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getLocalNodes (obj, num)
            if num > obj.TotalCells
                msg = ['\nCell number: ',num2str(num),' is greater than total number of cells: ',num2str(obj.TotalCells)];
                error('ErrorTests:convertTest',msg);
            end
            out = obj.ConnectivityArray{num};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getShiftedLocalNodes (obj, stride, num)
            n = obj.ConnectivityArray{num};
            out = uint32((stride(1))*ones(1,length(n))) + uint32(n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getFaceCellNodes (obj, fnum, cnum)
            out = obj.FaceCellNodes{fnum,cnum};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getShiftedFaceCellNodes (obj, stride, fnum, cnum)
            n = obj.FaceCellNodes{fnum,cnum};
            out = uint32(stride(1)*ones(1,length(n))) + uint32(n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate1DDoFs(obj, mesh)

verts = mesh.get_all_vertices();
if obj.FEMType == 1
    obj.TotalDoFs = obj.TotalCells * obj.Degree + 1;
elseif obj.FEMType == 2
    obj.TotalDoFs = obj.TotalCells * (obj.Degree + 1);
end
obj.NodeLocations = zeros(obj.TotalDoFs,obj.Dimension);

if obj.FEMType == 1
    for i=1:obj.TotalCells
        n = 1+(i-1)*(obj.Degree);
        obj.ConnectivityArray{i} = n:(n+obj.Degree);
        v = mesh.get_cell_verts(i);
        vv = sort(verts(v));
        dx = vv(2) - vv(1);
        ddx = dx/obj.Degree;
        obj.NodeLocations(n:(n+obj.Degree)) = vv(1):ddx:vv(2);
    end
elseif obj.FEMType == 2
    l = 0;
    for i=1:obj.TotalCells
        n = l+1;
        nn = l + obj.Degree + 1;
        obj.ConnectivityArray{i} = n:nn;
        v = mesh.get_cell_verts(i);
        vv = sort(verts(v));
        dx = vv(2) - vv(1);
        ddx = dx/obj.Degree;
        obj.NodeLocations(n:nn) = vv(1):ddx:vv(2);
        
        l = l + obj.Degree + 1;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate2DDoFs(obj, mesh)

verts = mesh.get_all_vertices();    nv = size(verts,1);
cells = mesh.get_all_cell_verts();
if obj.FEMType == 1
    if obj.Degree == 1
        obj.TotalDoFs = nv;
    elseif obj.Degree == 2
        if strcmp(obj.GeometryType,'Triangles') || strcmp(obj.GeometryType,'Quads')
            obj.TotalDoFs = nv + obj.TotalFaces;
        elseif strcmp(obj.GeometryType,'Quads')
            obj.TotalDoFs = nv + obj.TotalFaces + obj.TotalCells;
        elseif strcmp(obj.GeometryType,'Polygons')
            error('Cannot go higher than 1st order on polygons in 2D.')
        end
    else
        error('Cannot go higher than 2nd order in 2D.')
    end
elseif obj.FEMType == 2
    d = 0;
    if obj.Degree == 1
        for i=1:obj.TotalCells
            d = d + length(cells{i});
        end
    elseif obj.Degree == 2
        if strcmp(obj.GeometryType,'Triangles')
            for i=1:obj.TotalCells
                d = d + length(cells{i}) + length(mesh.get_cell_faces(i));
            end
        elseif strcmp(obj.GeometryType,'Quads')
            for i=1:obj.TotalCells
                d = d + length(cells{i}) + length(mesh.get_cell_faces(i)) + 1;
            end
        elseif strcmp(obj.GeometryType,'Polygons')
            error('Cannot go higher than 1st order on polygons in 2D.')
        end
    else
        error('Cannot go higher than 2nd order in 2D.')
    end
    obj.TotalDoFs = d;
end

obj.NodeLocations = zeros(obj.TotalDoFs,obj.Dimension);
% CFEM
if obj.FEMType == 1
    obj.NodeLocations(1:nv,:) = verts;
    obj.ConnectivityArray = cells;
    if obj.Degree == 2
        for i=1:obj.TotalCells
            n = length(obj.ConnectivityArray{i});
            f = mesh.get_cell_faces(i);
            nf = length(f);
            vv = nv*ones(1,nf) + f;
            obj.NodeLocations(vv,:) = mesh.get_face_centers(f);
            obj.ConnectivityArray{i}(n+1:n+nf) = vv;
            if strcmp(obj.GeometryType,'Quads')
                obj.ConnectivityArray{i}(n+nf+1) = nv+obj.TotalFaces+i;
                obj.NodeLocations(nv+obj.TotalFaces+i,:) = obj.get_cell_centers(i);
            end
        end
    end
    fcbool = ones(obj.TotalFaces,1);
    for i=1:obj.TotalCells
        f = mesh.get_cell_faces(i);
        nf = length(f);
        lca = obj.ConnectivityArray{i};
        cvs = mesh.get_cell_vert_indices(i);
        ncvs = size(cvs,1);
        for face=1:nf
            fvs = mesh.get_face_vert_indices(f(face));
            ffff = cvs(ismember(cvs,fvs));
            if obj.Degree == 2
                ffff = [ffff,lca(ncvs+face)];
            end
            obj.CellFaceNodes{i}{face} = ffff;
            obj.FaceCellNodes{f(face),fcbool(f(face))} = ffff;
            fcbool(f(face)) = fcbool(f(face)) + 1;
        end
    end
% DFEM
elseif obj.FEMType == 2
    % For the love of all that is holy,
    % please don't touch this...
    % Yes, I know it is not commented...
    nd = 0;
    fcbool = ones(obj.TotalFaces,1);
    for i=1:obj.TotalCells
        cf = mesh.get_cell_faces(i); nf = length(cf);
        cv = mesh.get_cell_vert_indices(i); nv = length(cv);
        ncd = nv;
        nodes = verts(cv,:);
        if obj.Degree == 2
            ncd = ncd + nf;
            nodes = [nodes;mesh.get_face_centers(cf)];
            if strcmp(obj.GeometryType,'Quads')
                ncd = ncd + 1;
                nodes = [nodes;mesh.get_cell_center(i)];
            end
        end
        for f=1:nf
            fverts = mesh.get_face_vert_indices(cf(f));
            for ff=1:length(fverts)
                ftemp = fverts(ff);
                fsame = find(cv==ftemp); % similar indices
                obj.FaceCellNodes{cf(f),fcbool(cf(f))} = [obj.FaceCellNodes{cf(f),fcbool(cf(f))},nd+fsame];
            end
            if obj.Degree == 2
                obj.FaceCellNodes{cf(f),fcbool(cf(f))} = [obj.FaceCellNodes{cf(f),fcbool(cf(f))},nd+nv+f];
            end
            obj.CellFaceNodes{i}{length(obj.CellFaceNodes{i})+1} = obj.FaceCellNodes{cf(f),fcbool(cf(f))};
            fcbool(cf(f)) = fcbool(cf(f)) + 1; % increment face index
        end
        obj.ConnectivityArray{i} = nd+1:nd+ncd;
        obj.NodeLocations(nd+1:nd+ncd,:) = nodes;
        nd = nd + ncd;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate3DDoFs(obj, mesh)

% Do this if the need arises, but it's going to suck...

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

