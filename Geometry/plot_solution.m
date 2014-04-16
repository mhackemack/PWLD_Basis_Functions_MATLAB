function plot_solution(mesh, DoF, x)

switch(mesh.Dimension)
    case(1)
        plot_1D_solution(mesh, DoF, x);
    case(2)
        plot_2D_solution(mesh, DoF, x);
    case(3)
        plot_3D_solution(mesh, DoF, x);
end