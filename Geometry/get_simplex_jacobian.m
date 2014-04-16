function J = get_simplex_jacobian(dim, verts)
J = zeros(dim,dim);
for i=1:dim
    J(:,i) = verts(i+1,:)' - verts(1,:)';
end