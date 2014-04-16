function plot_2D_mesh(mesh,coloring,PWLD)

hold on
% Loop through Cells
for e=1:mesh.TotalCells
    vverts = mesh.get_cell_verts(e);
    if nargin > 1
        if coloring == 1
            fill(vverts(:,1),vverts(:,2),get_color_map(mesh.get_cell_material_id(e)));
        else
            fill(vverts(:,1),vverts(:,2),[1 1 1]);
        end
    else
        fill(vverts(:,1),vverts(:,2),[1 1 1]);
    end
    if nargin == 3
        if logical(PWLD)
            rcenter = mean(vverts);
            plot(rcenter(1),rcenter(2),'ok')
            for i=1:size(vverts,1)
                plot([rcenter(1),vverts(i,1)],[rcenter(2),vverts(i,2)],'--k')
            end
        end
    end
end
% Loop through Edges
for e=1:mesh.TotalFaces
    vverts = mesh.get_face_verts(e);
    if nargin == 3
        if logical(PWLD)
            plot(vverts(:,1),vverts(:,2),'k','LineWidth',2);
        else
            plot(vverts(:,1),vverts(:,2),'k');
        end
    else
        plot(vverts(:,1),vverts(:,2),'k');
    end
end
hold off

function out = get_color_map(val)

switch(val)
    case (1)
        out = [0 0 1];
    case (2)
        out = [0 1 0];
    case (3)
        out = [1 0 0];
    case (4)
        out = [1 0 1];
    case (5)
        out = [0 1 1];
    case (6)
        out = [1 1 0];
    otherwise
        out = [0 0 1];
end
