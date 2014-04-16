%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Discontinuous (PWLD) Basis Function 
%                   Generator - Surface Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's surface using 
%                   the PWLD DGFEM basis functions.
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
function [varargout] = PWLD_surface(verts, faces)

if nargin == 0
    error('--- No inputs specified. ---')
else
    % Prepare Vertices and Dimensional Space
    % --------------------------------------
    [mv,nv] = size(verts);
    if nv > mv, verts = verts'; end
    [nv,dim] = size(verts);
    rcenter = get_center_point(verts);
    
    % Allocate Matrix Memory
    % ----------------------
    M = zeros(nv,nv);           % mass matrix
    G = cell(dim,1);            % gradient matrix
    for i=1:dim
        G{i} = zeros(nv,nv);
    end
    
    % 1D Generation
    % -------------
    if dim == 1
        error('Choosing not to do PWLD in 1D -- is this just LD???')
    end
    % 2D Generation
    % -------------
    if dim == 2
        ffaces{1} = 1:nv;
    end
    % 3D Generation
    % -------------
    if dim == 3
        % Check face structure
        % --------------------
        if nargin ~= 2
            error('--- No face input specified. ---')
        end
        if isempty(faces)
            error('--- No face input specified. ---')
        end
        if iscell(faces)
            ffaces = faces;
        else
            ffaces = cell(size(faces,1),1);
            for i=1:length(ffaces)
                ffaces{i} = faces(i,:);
            end
        end
    end
    % Loop through Faces
    % ------------------
    for f=1:length(ffaces)
        ff = ffaces{f};
        fcenter = get_center_point(verts(ff,:));
        ne = length(ff);
        % Loop through Edges on Face
        % --------------------------
        for e=1:ne
            if e == ne
                ee = ff(1);
            else
                ee = ff(e+1);
            end
            eee = [ff(e),ee];
            % Edge Triangle/Tetrahedron Information
            % -------------------------------------
            if dim==2
                lverts = [verts(eee,:);rcenter];
            else
                lverts = [verts(eee,:);fcenter;rcenter];
            end
            [J,invJ,detJ,Vside] = get_jacobian(lverts);
            [lens,vecs] = get_side_lengths(lverts);
            % Reference Triangle/Tetrahedron Matrices
            % ---------------------------------------
            m = get_ref_mass_matrix(dim)*Vside;
%             s = get_local_stiffness_matrix(dim,lens,Vside);
            g = get_local_gradient_term(dim,lens,vecs,Vside,invJ,detJ);
            % Append to Global Matrices
            % -------------------------
            M = M + matrix_contribution(dim,nv,m,eee,ff);
%             K = K + matrix_contribution(dim,nv,s,eee,ff);
            for j=1:dim
                G{j} = G{j} + matrix_contribution(dim,nv,g{j},eee,ff);
            end
        end
    end
    % Set Outputs
    varargout{1} = M;
    varargout{2} = G;
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_center_point(verts)
out = mean(verts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, invJ, detJ, vol] = get_jacobian(verts)
dim = size(verts,2);
J = zeros(dim,dim);
for i=1:dim
    J(:,i) = verts(i+1,:)' - verts(1,:)';
end
invJ = J^(-1);
detJ = abs(det(J));
vol = detJ/factorial(dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grads(dim)
out = [-ones(1,dim);diag(ones(dim,1))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_mass_matrix(dim)
if dim==2
    out = [2,1,0;1,2,0;0,0,0]./6;
elseif dim==3
    out = [2,1,1,0;1,2,1,0;1,1,2,0;0,0,0,0]./12;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_local_stiffness_matrix(dim,lens,A)
R = lens.^2/(4*A);
if dim == 2
    
elseif dim == 3
    out = [  2*R(1), R(3) - R(1) - R(2), R(2) - R(1) - R(3);...
             R(3) - R(1) - R(2), 2*R(2), R(1) - R(2) - R(3);...
             R(2) - R(1) - R(3), R(1) - R(2) - R(3), 2*R(3)    ]./2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_local_gradient_term(dim,lens,vecs,vol,invJ,detJ)
out = cell(dim,1);
if dim == 2
    a = -lens(end)*([lens,lens].*vecs)./(2*vol);
    b = [1/2,0,0];
elseif dim == 3
%     a = -lens(end)*([lens,lens,lens].*vecs)./(6*vol);
    b = [1/6,1/6,0,0];
    db = get_basis_grads(dim);
    a = detJ*db*invJ;
end
for i=1:dim
    out{i} = a(:,i)*b;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lens,vecs] = get_side_lengths(verts)
[nv,dim] = size(verts);
lens = zeros(nv,1);
vecs = zeros(nv,dim);
if dim == 2
    dd = verts(3,:)-verts(2,:); lens(1) = norm(dd); vecs(1,:) = [dd(2),-dd(1)]./lens(1);
    dd = verts(1,:)-verts(3,:); lens(2) = norm(dd); vecs(2,:) = [dd(2),-dd(1)]./lens(2);
    dd = verts(2,:)-verts(1,:); lens(3) = norm(dd); vecs(3,:) = [dd(2),-dd(1)]./lens(3);
elseif dim == 3
    ord = get_tet_ordering();
    for i=1:4
        tord = ord(i,:);
        u = cross([verts(tord(1),:)-verts(tord(2),:)]',[verts(tord(3),:)-verts(tord(2),:)]')';
        lens(i) = norm(u,2)/2;
        vecs(i,:) = u./norm(u);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = matrix_contribution(dim,nv,mat,v,fv)
a = 1/nv;
out = zeros(nv,nv);
out(v,v) = mat(1:length(v),1:length(v));
for i=1:length(v)
    out(v(i),:) = out(v(i),:) + a*mat(i,end);
    out(:,v(i)) = out(:,v(i)) + a*mat(end,i);
end
out(:,:) = out(:,:) + a*a*mat(end,end);
if dim == 3
    b = 1/length(fv);
    out(fv,:) = out(fv,:) + a*b*mat(end-1,end);
    out(:,fv) = out(:,fv) + a*b*mat(end,end-1);
    out(fv,fv) = out(fv,fv) + b*b*mat(end-1,end-1);
    for i=1:length(v)
        out(v(i),fv) = out(v(i),fv) + b*mat(i,end-1);
        out(fv,v(i)) = out(fv,v(i)) + b*mat(end-1,i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_tet_ordering()
out = [2,3,4;1,4,3;1,2,4;1,3,2];