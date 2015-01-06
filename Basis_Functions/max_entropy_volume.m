%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy Basis Function Generator - Volume Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Maximum-Entropy basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes - 2D:     1) 'faces' input not needed (ignored)
%                   2) 'verts' input is in the form (npts x ndim)
%                   3) 'verts' coordinates need to be CCW
%
%   Notes - 3D:     1) 'verts' and 'faces' input both needed
%                   2) 'faces' holds the vertex numberings of 'verts'
%                   3) 'faces' can take either cell or array structure 
%                      - if array structure: form (nfaces x npts_face)
%                                            where npts_face is constant
%                      - array structure can only be used if each face
%                        has the same number of vertices
%                   4) Vertices in 'verts' do not need any proper ordering
%                   5) Vertices on each face need to be in CCW order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = max_entropy_volume(varargin)
nout = nargout;
% Collect Input Arguments
% -----------------------
nverts = varargin{5};
verts = varargin{1}(1:nverts,:);
faces = varargin{2};
flags = varargin{3};
order = varargin{4};
% Quick Error Checking
% --------------------
if order > 2, error('2 is the maximum order (serendipity).'); end
% Prepare Vertices and Dimensional Space
% --------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
dim = size(verts, 2);
% Allocate Matrix Memory
% ----------------------
if dim == 2
    ntot = order*nverts;
    nfaces = nverts;
else
    nfaces = length(faces);
    ntot = nverts + (order-1)*(nverts+nfaces-2); % Euler's formula
end
M = zeros(ntot,ntot);       % mass matrix
K = zeros(ntot,ntot);       % stiffness matrix
G = cell(dim,1);            % gradient matrix
for i=1:dim
    G{i} = zeros(ntot,ntot);
end
% Get Problem Preliminaries
% -------------------------
[qx, qw] = get_general_volume_quadrature(verts, faces, 2*order+1); nqx = length(qw);
rva = get_vertex_differences(verts, qx);

% Newton Iterations
% -----------------
iter = 1; converged = false; 
x = zeros(dim*nqx, 1);
g = get_F( x, order );
while ~converged
    % Get Hessian
    H = get_J( x, order );
    % Get search direction
    dx = -H\g;
    % Get step length
    a = get_dampen_term( dx, order );
    % Update lambda multipliers
    x = x + a*dx;
    % Check convergence
    g = get_F( x );
    if norm(g) < 1e-12; converged = true; end
end

% Perform Matrix Generations
% --------------------------
Z = zeros(nqx, ntot);
% Compute Gradients
if flags == 2 || flags == 3
    
end
% Switch method order


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nqx = size(qx, 1);
out = zeros(nv, dim, nqx);
for q=1:nqx
    for d=1:dim
        out(:,d,q) = verts(:,d) - qx(q,d);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_F( v, order )

if order == 1
    
elseif order == 2
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_J( v, order )

if order == 1
    
elseif order == 2
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_dampen_term( v, order )
out = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
