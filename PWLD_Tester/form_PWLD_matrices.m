function [M,K,g] = form_PWLD_matrices(verts, faces)

[M,K,~] = PWLD_volume(verts, faces);
[~,~,g] = PWLD_surface(verts, faces);

