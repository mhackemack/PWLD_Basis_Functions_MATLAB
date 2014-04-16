function [x,w,opp_ind] = get_angular_quadrature(dim, AQName, level, npolar, nazimuth)

if dim == 1
    [x,w] = lgwt(level,-1,1);
    w = w*2;
else
    if strcmp(AQName, 'LDFE')
        [x,w] = get_LDFE_quad(dim, level);
    elseif strcmp(AQName, 'PGLC')
        [x,w] = get_PGLC_quad(dim, npolar, nazimuth);
    elseif strcmp(AQName, 'LS')
        [x,w] = get_LS_quad(dim, level);
    elseif strcmp(AQName, 'TriGLC')
        [x,w] = get_TriGLC_quad(dim, npolar, nazimuth);
    end
end

[x,w,opp_ind] = deploy_all_octants(dim,x,w);
x = x(:,1:dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_LDFE_quad(dim, level)
octAngles = 4^(level+1);
numAngles = octAngles*2^dim;


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_PGLC_quad(dim, npolar, nazimuth)
numAngles = npolar*nazimuth*2^dim;
nRoots = 2*npolar;
[x,w] = lgwt(nRoots,-1,1);

% Allocate output memory
angs = zeros(numAngles,3);
wts = zeros(numAngles,1);

% Get other local variables
delta_phi = pi/(4*nazimuth);
azim_weight = 2*delta_phi;
azim_angle = delta_phi;

for i=1:nazimuth
    for j=npolar-1:nRoots
        costheta = x(j);
        sintheta = sqrt(1-costheta^2);
        iord = j + (i-1)*npolar;
        angs(iord,:) = [cos(azim_angle)*sintheta,sin(azim_angle)*sintheta,costheta];
        wts(iord) = azim_weight*w(j);
    end
    azim_angle = azim_angle + 2*delta_phi;
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_LS_quad(dim, level)
if level > 24
    error('Error: Cannot go high thatn S24 for level-symmetric. Yields negative weights.')
end
[angs,wts] = get_LS_local_quad(dim,level);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_TriGLC_quad(dim, npolar, nazimuth)
octAngles = 0;
nnazimuth = nazimuth;
for i=1:npolar
    octAngles = octAngles + nnazimuth;
    nnazimuth = nnazimuth - 1;
end
numAngles = octAngles*2^dim;
nRoots = 2*npolar;
[x,w] = lgwt(nRoots,-1,1);
% Allocate memory
angs = zeros(numAngles,3);
wts  = zeros(numAngles,1);
% Loop through polar directions
dir = 1;
for j=1:npolar
    nazimuth_ = nazimuth;
    delta_phi = pi/(4.0*nazimuth_);
    azim_angle = delta_phi;
    azim_weight = 2*delta_phi;
    
    k = j + npolar;
    costheta = x(k);
    sintheta = sqrt(1-costheta^2);
    for i=1:nazimuth
        angs(dir,:) = [cos(azim_angle)*sintheta,sin(azim_angle)*sintheta,costheta];
        wts(dir) = azim_weight*w(k);
        azim_angle = azim_angle + 2*delta_phi;
        dir = dir + 1;
    end
    nazimuth = nazimuth - 1;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts,opp_ind] = deploy_all_octants(dim, angs, wts)
numAngles = length(wts);
octAngles = numAngles / 2^dim;
opp_ind = zeros(numAngles,1);
if dim == 1
    for n=1:octAngles
        opp_ind(n) = n + 2*octAngles;
        opp_ind(n + 2*octAngles) = n;
    end
end
if dim == 2 || dim == 3
    for octant=2:4
        for n=1:octAngles
            m = (octant - 1)*octAngles + n;
            switch (octant)
                case(2)
                    angs(m,1) = -angs(n,1);
                    angs(m,2) =  angs(n,2);
                    angs(m,3) =  angs(n,3);
                case(3)
                    angs(m,1) = -angs(n,1);
                    angs(m,2) = -angs(n,2);
                    angs(m,3) =  angs(n,3);
                case(4)
                    angs(m,1) =  angs(n,1);
                    angs(m,2) = -angs(n,2);
                    angs(m,3) =  angs(n,3);
            end
            wts(m) = wts(n);
        end
    end
end
if dim == 3
    for n=1:numAngles/2
        angs(n+numAngles/2,1) =  angs(n,1);
        angs(n+numAngles/2,2) =  angs(n,2);
        angs(n+numAngles/2,3) = -angs(n,3);
        wts(n+numAngles/2) = wts(n);
    end
end
if dim == 2
    for n=1:octAngles
        opp_ind(n) = n+2*octAngles;
        opp_ind(n+2*octAngles) = n;
        opp_ind(n+octAngles) = n+3*octAngles;
        opp_ind(n+3*octAngles) = n+octAngles;
    end
elseif dim == 3
    for n=1:octAngles
        opp_ind(n) = n+6*octAngles;
        opp_ind(n+6*octAngles) = n;
        opp_ind(n+octAngles) = n+7*octAngles;
        opp_ind(n+7*octAngles) = n+octAngles;
        opp_ind(n+2*octAngles) = n+4*octAngles;
        opp_ind(n+4*octAngles) = n+2*octAngles;
        opp_ind(n+3*octAngles) = n+5*octAngles;
        opp_ind(n+5*octAngles) = n+3*octAngles;
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

