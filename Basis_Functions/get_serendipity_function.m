function out = get_serendipity_function(flag, ctype, dim)

if strcmp(flag, 'vals')
    if dim == 1
        out = @evaluate_1D_reference_values;
    elseif dim == 2
        if ctype == 1
            out = @evaluate_2D_triangle_reference_values;
        else
            out = @evaluate_2D_quad_reference_values;
        end
    else
        if ctype == 1
            out = @evaluate_3D_tet_reference_values;
        else
            out = @evaluate_3D_hex_reference_values;
        end
    end
elseif strcmp(flag, 'grads')
    if dim == 1
        out = @evaluate_1D_reference_gradients;
    elseif dim == 2
        if ctype == 1
            out = @evaluate_2D_triangle_reference_gradients;
        else
            out = @evaluate_2D_quad_reference_gradients;
        end
    else
        if ctype == 1
            out = @evaluate_3D_tet_reference_gradients;
        else
            out = @evaluate_3D_hex_reference_gradients;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Basis Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_values(deg, x)
if deg == 1
    out = [1-x, x];
elseif deg == 2
    out = [2*(1/2-x).*(1-x), -2*x.*(1/2-x), 4*x.*(1-x)];
elseif deg == 3
    out = [9/2*(1/3-x).*(2/3-x).*(1-x), 9/2*(1/3-x).*(2/3-x).*x, 27/2*(1-x).*(2/3-x).*x, -27/2*(1-x).*(1/3-x).*x];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_gradients(deg, x)
if deg == 1
    out = [-1, 1];
elseif deg == 2
    out = [4*x-3, 4*x-1, 4-8*x];
elseif deg == 3
    out = [-9/2*(3*x.^2 - 4*x + 11/9), 9/2*(3*x.^2 - 2*x + 2/9), 9/2*(9*x.^2 - 10*x + 2), -9/2*(9*x.^2 - 8*x + 1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_values(deg, x)
s = x(:,1); t = x(:,2);
if deg == 1
    out = [1-s-t,s,t];
elseif deg == 2
    out = [  2*s.^2+4*s.*t+2*t.^2-3*s-3*t+1,...
             2*s.^2-s,...
             2*t.^2-t,...
             -4*s.^2-4*s.*t+4*s,...
             4*s.*t,...
             -4*t.^2-4*s.*t+4*t...
          ];
elseif deg == 3
    out = [(1-s-t).*(1-3*s-3*t).*(1-3/2*s-3/2*t),...
            s.*(1-3*s).*(1-3/2*s),...
            t.*(1-3*t).*(1-3/2*t),...
            9*s.*(1-s-t).*(1-3/2*s-3/2*t),...
            -9/2*s.*(1-s-t).*(1-3*s),...
            -9/2*s.*t.*(1-3*s),...
            -9/2*s.*t.*(1-3*t),...
            -9/2*(1-s-t).*t.*(1-3*s),...
            9*(1-s-t).*t.*(1-3/2*s-3/2*t)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_values(deg, xx)
s = xx(:,1); t = xx(:,2);
if deg == 1
    out = [(1-s).*(1-t), s.*(1-t), s.*t, t.*(1-s)];
elseif deg == 2
    out = [(1-s).*(1-t).*(1-2*s-2*t),...
        -s.*(1-t).*(1-2*s+2*t),...
        -s.*t.*(3-2*s-2*t),...
        -(1-s).*t.*(1+2*s-2*t),...
        4*s.*(1-s).*(1-t),...
        4*s.*t.*(1-t),...
        4*s.*(1-s).*t,...
        4*t.*(1-s).*(1-t)];
elseif deg == 3
    out = []./2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_gradients(deg, xx)
if deg == 1
    out = [-1,-1;1,0;0,1];
elseif deg == 2
    s = xx(:,1); t = xx(:,2);
    out = [  4*s+4*t-3,  4*s+4*t-3;...
             4*s-1,         0;...
             0,                4*t-1;...
             -8*s-4*t+4,-4*s;...
             4*t,           4*s;...
             -4*t,    -8*t-4*s+4];
elseif deg == 3
    out = [1/2-3/2*(3*s+3*t-2).^2,1/2-3/2*(3*s+3*t-2).^2;...
           1+9/2*s.*(3*s-2),0;...
           0,1+9/2*t.*(3*t-2);...
           81/2*s.*s+54*s.*t-45*s+27/2*t.*t-45/2*t+9,27*s.*t+27*s.*s-45/2*s;...
           -81/2*s.*s-9/2*s.*(6*t-8)+9/2*t-9/2,9/2*s-27/2*s.*s;...
           27*s.*t-9/2*t,27/2*s.*s-9/2*s;...
           27*t.*t-9/2*t,27*s.*t-9/2*s;...
           9/2*t-27/2*t.*t,-27*s.*t+9/2*s-81/2*t.*t+36*t-9/2;...
           27*s.*t+27*t.*t-45/2*t,81/2*t.*t+54*s.*t-45*t+27/2*s.*s-45/2*s+9];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_gradients(deg, xx)
s = xx(:,1); t = xx(:,2);
if deg == 1
    out = [t-1,s-1;1-t,-s;t,s;-t,1-s];
elseif deg == 2
    out = [-3+5*t-4*s.*t-2*t.*t+4*s,-3+5*s-4*s.*t-2*s.*s+4*t;...
           4*s-4*s.*t+2*t.*t-t-1,-2*s.*s+4*s.*t-s;...
           4*s.*t+2*t.*t-3*t,2*s.*s+4*s.*t-3*s;...
           4*s.*t-2*t.*t-t,2*s.*s-4*s.*t-s-1+4*t;...
           8*s.*t-8*s-4*t+4,4*s.*s-4*s;...
           4*t-4*t.*t,4*s-8*s.*t;...
           4*t-8*s.*t,4*s-4*s.*s;...
           4*t.*t-4*t,8*s.*t-4*s-8*t+4];
elseif deg == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_tet_reference_values(deg, xx)
x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    
elseif deg == 2
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_tet_reference_gradients(deg, xx)
x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    
elseif deg == 2
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_hex_reference_values(deg, xx)
x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    out = [(1-x).*(1-y).*(1-z),...
            x.*(1-y).*(1-z),...
            x.*y.*(1-z),...
            (1-x).*y.*(1-z),...
            (1-x).*(1-y).*z,...
            x.*(1-y).*z,...
            x.*y.*z,...
            (1-x).*y.*z];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_hex_reference_gradients(deg, xx)
x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    out = [-y.*z+y+z-1,  -x.*z+x+z-1,  -x.*y+x+y-1;...
           (1-y).*(1-z),  x.*z-x,       x.*y-x;...
           y.*(1-z),      x.*(1-z),    -x.*y;...
           y.*z-y,       (1-x).*(1-z),  x.*y-y;...
           y.*z-z,       x.*z-z         (1-x).*(1-y);...
           (1-y).*z,     -x.*z          x.*(1-y);...
           y.*z,          x.*z,         x.*y;...
           -y.*z,         (1-x).*z,     (1-x).*y];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%