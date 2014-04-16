classdef FEHandler < handle
    properties (Access = public)
        Dimension
        Degree
        FEMType
        FEMName
    end
    properties (Access = private)
        TotalDoFs
        TotalCells
        TotalFaces
        
        CurrentCell
        GlobalShift
    end
    properties (Access = private)
        CellVertices
        Jacobians
        InverseJacobians
        StartingPositions
    end
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                           Constructor
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = FEHandler (varargin)
            n = nargin;
            if n == 0
                % empty constructor -> do nothing
            elseif n == 1
                error('Not enough input arguments...')
            else
                if isa(varargin{1}, 'DoFHandler')
                    mesh = varargin{2};
                    DoF = varargin{1};
                else
                    mesh = varargin{1};
                    DoF = varargin{2};
                end
                clear varargin
                
                obj.Dimension = mesh.Dimension;
                obj.TotalCells = mesh.TotalCells;
                obj.TotalFaces = mesh.TotalFaces;
                obj.Degree = DoF.Degree;
                obj.TotalDoFs = DoF.TotalDoFs;
                
                % Generate Memory Space
                obj.Jacobians = cell(obj.TotalCells,1);
                obj.InverseJacobians = cell(obj.TotalCells,1);
                obj.StartingPositions = zeros(obj.TotalCells,obj.Dimension);
                obj.CellVertices = cell(obj.TotalCells,1);
                
                for c=1:obj.TotalCells
                    obj.CellVertices{c} = mesh.get_cell_verts(c);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update_cell(obj, cellID)
            obj.CurrentCell = cellID;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = private)
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_values(obj, x)
if obj.Degree == 1
    out = [(1-x)/2,(x+1)/2];
elseif obj.Degree == 2
    out = [x*(x-1)/2,1-x*x,x*(x+1)/2];
elseif obj.Degree == 3
    out = [  9*(1/9-x*x)*(x-1)/16, ...
             27*(1-x*x)*(1/3-x)/16,...
             27*(1-x*x)*(1/3+x)/16,...
             -9*(1/9-x*x)*(1+x)/16 ...
             ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_gradients(obj, x)
if obj.Degree == 1
    out = [-1/2;1/2];
elseif obj.Degree == 2
    out = [x-1/2;-2*x;x+1/2];
elseif obj.Degree == 3
    out = [  -9*(3*x*x-2*x-1/9)/16;...
             27*(3*x*x-2*x/3-1)/16;...
             27*(-3*x*x-2*x/3+1)/16;...
             -9*(-3*x*x-2*x+1/9)/16 ...
             ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_values(obj, x)
if obj.Degree == 1
    out = [1-x(1)-x(2),x(1),x(2)];
elseif obj.Degree == 2
    out = [  2*x(1)^2+4*x(1)*x(2)+2*x(2)^2-3*x(1)-3*x(2)+1,...
             2*x(1)^2-x(1),...
             2*x(2)^2-x(2),...
             -4*x(1)^2-4*x(1)*x(2)+4*x(1),...
             4*x(1)*x(2),...
             -4*x(1)^2-4*x(1)*x(2)+4*x(1)...
          ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_values(obj, x)
if obj.Degree == 1
    out = 1/4*[(1-x(1))*(1-x(2)),(1+x(1))*(1-x(2)),(1+x(1))*(1+x(2)),(1-x(1))*(1+x(2))];
elseif obj.Degree == 2
    out = [  (1-x(1))*(x(2)-1)*(x(1)+x(2)+1)/4,...
             (1+x(1))*(x(2)-1)*(x(1)+x(2)+1)/4,...
             (1+x(1))*(x(2)+1)*(x(1)+x(2)-1)/4,...
             (x(1)-1)*(x(2)+1)*(x(1)-x(2)+1)/4,...
             (1-x(2))*(1-x(1)^2)/2,...
             (1+x(1))*(1-x(2)^2)/2,...
             (1+x(2))*(1-x(1)^2)/2,...
             (1-x(1))*(1-x(2)^2)/2                 ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_gradients(obj, x)
if obj.Degree == 1
    out = [-1,-1;1,0;0,1];
elseif obj.Degree == 2
    out = [  4*x(1)+4*x(2)-3,  4*x(1)+4*x(2)-3;...
             4*x(1)-1,         0;...
             0,                4*x(2)-1;...
             -8*x(1)-4*x(2)+4,-4*x(1);...
             4*x(1),           4*x(2);...
             -8*x(1)-4*x(2)+4,-4*x(1)    ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_gradients(obj, x)
if obj.Degree == 1
    out = 1/4*[-1+x(2), -1+x(1);...
                1-x(2), -1-x(1);...
                1+x(2),  1+x(1);...
               -1-x(2),  1-x(1)];
elseif obj.Degree == 2
    out = [  -1/4*(x(2)-1)*(2*x(1)+x(2)),   -1/4*(x(1)-1)*(x(1)+2*x(2));...
             -1/4*(x(2)-1)*(2*x(1)-x(2)),   -1/4*(x(1)+1)*(x(1)-2*x(2));...
              1/4*(x(2)+1)*(2*x(1)+x(2)),    1/4*(x(1)+1)*(x(1)+2*x(2));...
              1/4*(x(2)+1)*(2*x(1)-x(2)),    1/4*(x(1)-1)*(x(1)-2*x(2));...
              x(1)*(x(2)-1),                 1/2*(x(1)^2-1);...
              1/2*(1-x(2)^2),                -(x(1)+1)*x(2);...
              -(x(2)+1)*x(1),                1/2*(1-x(1)^2);...
              1/2*(x(2)^2-1),                (x(1)-1)*x(2)...
          ];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = evaluate_2D_PWLD_value(obj, x)
% cv = obj.CellVertices{obj.CurrentCell};
% ccent = mean(cv);
% nv = size(cv,1); a = 1/nv;
% for v=1:nv
%     if v==nv
%         vv = [v,1];
%     else
%         vv = [v,v+1];
%     end
%     vvv = [cv(vv);ccent];
%     if logical(inpoly(x,vvv))
%         J = get_simplex_jacobian(obj.Dimension,vvv);
%     end
% end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = evaluate_3D_PWLD_value(obj, x)
% cv = obj.CellVertices{obj.CurrentCell};
% ccent = mean(cv);
% nv = size(cv,1); a = 1/nv;
% 
% 
% end