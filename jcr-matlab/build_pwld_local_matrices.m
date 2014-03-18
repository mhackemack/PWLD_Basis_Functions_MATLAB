function [M,K,f,g]=build_pwld_local_matrices(poly,vert)
% vertices are entered anti-clockwise
% poly=[1 2 3]
% vert=[0 0; 1 0; 0 1]

global verbose

% centroid
C=mean(vert);
% alpha coef
nv=length(poly);
alpha=1/nv;

% mass and stiffness matrices for the entire polygon
M=zeros(nv,nv);   % mass mat
K=zeros(nv,nv);   % stiffness mat
f=zeros(nv,1);    % rhs vector
% mass integrals over one side:
M_side=zeros(nv,nv);
K_side=zeros(nv,nv);
f_side=zeros(nv,1);
% gradient on all sides
g=zeros(2,nv,nv); % x/y | tf | iside
% integrals of b^i_1/2 b^i_1/2 :
m=[2 1;1 2]/24+alpha*(1+alpha)*[1 1;1 1]/12;
% integrals of b^i_1/2 b^i_c :
aux1=alpha/24*(1+2*alpha); % a/24 + a^2/12
% integrals of b^i_c b^i_c :
aux2=alpha^2/12;
% mass matrix
M_side([1 2],[1 2]) = m;
M_side([1 2],3:end) = aux1;
M_side(3:end,[1 2]) = aux1;
M_side(3:end,3:end) = aux2;
% rhs
f_side(1:2)=(1+alpha)/6;
f_side(3:end)=alpha/6;
% gradient per side
g_side=[-1 1 0; (alpha-1) alpha alpha];
% list of vertices, when looping over sides
list_vert=1:nv;

% loop over sides
for iside=1:nv
    
    % pick 1st vertex
    irow1=iside;
    % pick 2ndt vertex
    irow2=irow1+1; if(irow2>nv), irow2=1; end
    
    % compute |J|=2.Area for that side
    A=vert(irow1,:); B=vert(irow2,:);
    Jac_i = cross([B-A 0]',[C-A 0]'); % AB ^ AC
    if verbose
        if(Jac_i(3)<0),
            [A B C]
            poly
            vert
            iside
            warning('negative jac, ordering of the side is not clockwise, indicative of centroid being outside of polygonal the mesh');
        end
    end
    det_J_i=norm(Jac_i,2);
    %%% hack
    if Jac_i(3)<0, det_J_i=-det_J_i; end
    %%% end hack
    Jac_i=[(B-A)' (C-A)'];
%     det_J_i=det(Jac_i); % safer
    iJt=inv(Jac_i');
    if verbose
        if(abs(det(Jac_i)-det_J_i)>1e-11),
            abs(det(Jac_i)-det_J_i)
            fprintf('Jac1= %E /= Jac2= %E \n',det(Jac_i),det_J_i);
            warning('2 versions of Jac do not yield same det ...');
        end
    end
    % add contribution from M_side
    M(list_vert,list_vert) = M(list_vert,list_vert) + det_J_i*M_side;
    
    % add contribution from f_side
    f(list_vert) = f(list_vert) + det_J_i*f_side;
    
    % stiffness matrix
    a = alpha;
    L1=norm(B-A); L2=norm(C-B); L3=norm(A-C);
    r1 = L1^2/(2*det_J_i);
    r2 = L2^2/(2*det_J_i);
    r3 = L3^2/(2*det_J_i);
    % [r1 r2 r3]
    % sanity check:
    test = 2*r1*(r2+r3)+2*r2*r3-(r1^2+r2^2+r3^2)-1;
    if verbose
        if(abs(test)>1e-10), test, warning('r_i error in stiffness matrix coefficients'); end
    end
    % stiffness matrix for side i
    kk_side=[ ...
        (-1 + a)*a*r1 - (-1 + a)*r2 + a*r3,  ((1 - 2*a + 2*a^2)*r1 - r2 - r3)/2,  (a*((-1 + 2*a)*r1 - r2 + r3))/2;
        ((1 - 2*a + 2*a^2)*r1 - r2 - r3)/2,  (-1 + a)*a*r1 + a*r2 - (-1 + a)*r3,  (a*((-1 + 2*a)*r1 + r2 - r3))/2;
        (a*((-1 + 2*a)*r1 - r2 + r3))/2   ,  (a*((-1 + 2*a)*r1 + r2 - r3))/2,     a^2*r1
        ];
    K_side([1 2],[1 2]) = kk_side([1 2],[1 2]);
%     K_side([1 2],3:end) = kk_side(1,3);
%     K_side(3:end,[1 2]) = kk_side(1,3);
    K_side(1,3:end) = kk_side(1,3);
    K_side(2,3:end) = kk_side(2,3);
    K_side(3:end,1) = kk_side(3,1);
    K_side(3:end,2) = kk_side(3,2);
    K_side(3:end,3:end) = kk_side(3,3);
    % add contribution from K_side
    K(list_vert,list_vert) = K(list_vert,list_vert) + K_side;
    
    % gradient 
    aux=iJt*g_side;
    g(:,list_vert(1:2),iside) = aux(1:2,1:2);
    for k=3:nv
        g(:,list_vert(k),iside) = aux(1:2,3);
    end
    
    
    % shift vertex IDs
    list_vert=[list_vert list_vert(1)];
    list_vert(1)=[];
end

return
end
