function [sol, geometry, DoF] = execute_problem(data, geometry)

% Generate problem execution schedule
% -----------------------------------
if strcmp(data.problem.problemType,'SourceDriven')
    pcall = @perform_source_driven_problem;
elseif strcmp(data.problem.problemType,'Eigenvalue')
    pcall = @perform_eignevalue_problem;
end
if data.problem.refineMesh == 0
    piters = 1;
    rfunc = @emptyFunction;
else
    if data.problem.refinementLevels < 1
        piters = 1;
        rfunc = @emptyFunction;
    else
        piters = data.problem.refinementLevels;
        rfunc = @refine_mesh;
    end
end
% Execute problem based on user input
% -----------------------------------
for iter=1:piters
    DoF = DoFHandler(geometry, data.Neutronics.FEMDegree, data.Neutronics.FEMType);
    sol = solution_allocation(data.Neutronics, DoF);
    sol = pcall(data, geometry, DoF, sol);
    iters{iter} = sol.iter;
    times{iter} = sol.times;
    % Refine mesh if necessary
    rfunc(geometry, DoF, data.problem, sol.flux);
end

% Set problem outputs
sol.iters = iters;
sol.times = times;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = perform_source_driven_problem(data, geometry, DoF, sol)

% Create DoF and other structures
% -------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
fcall = get_solution_function_handle(lower(data.Neutronics.FEMType), data.Neutronics.transportMethod);
sol.flux0 = sol.flux;

% Loop through iterations
% -----------------------
for l=1:data.solver.maxIterations
    disp(['Perform Source-Driven Iteration: ',num2str(l)])
    
    tictime = tic;
    sol.flux = fcall(ndat,solvdat,geometry,DoF,sol.flux);
    ferr = compute_flux_moment_differences(sol.flux,sol.flux0,1:sol.numberEnergyGroups,1);
    % Output iteration data
    disp(['Flux Error: ',num2str(ferr)])
    disp(' ')
    % Check for Convergence
    if ferr < solvdat.relativeTolerance
        toctime(l) = toc(tictime);
        break
    else
        toctime(l) = toc(tictime);
        sol.flux0 = sol.flux;
    end
end

% Apply outputs
% -------------
sol.iter = l;
sol.times = toctime;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = perform_eignevalue_problem(data, geometry, DoF, sol)

% Create DoF and other structures
% -------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
fcall = get_solution_function_handle(lower(data.Neutronics.FEMType), data.Neutronics.transportMethod);

% Loop through iterations
% -----------------------
keff0 = ndat.keff;
sol.flux0 = sol.flux;
for l=1:data.solver.maxIterations
    disp(['Perform Eigenvalue Iteration: ',num2str(l)])
    
    tictime = tic;
    sol.flux = fcall(ndat,solvdat,geometry,DoF,sol.flux);
    % Compute new keff and errors
    [keff,sol.flux] = estimate_new_keff(sol.flux,sol.flux0);
    ferr = compute_flux_moment_differences(sol.flux,sol.flux0,1:sol.numberEnergyGroups,1);
    kerr = abs(keff - keff0) / abs(keff);
    % Output iteration data
    disp(['keff: ',num2str(keff)])
    disp(['Flux Error: ',num2str(ferr)])
    disp(['keff Error: ',num2str(kerr)])
    disp(' ')
    % Check for Convergence
    if ferr < solvdat.relativeTolerance && kerr < solvdat.relativeTolerance
        toctime(l) = toc(tictime);
        break
    else
        toctime(l) = toc(tictime);
        sol.flux0 = sol.flux;
        keff0 = keff;
        ndat.keff = keff;
    end
    
end

% Apply outputs
% -------------
sol.keff = keff;
sol.iter = l;
sol.times = toctime;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      Miscellaneous Function Calls
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcall = get_solution_function_handle(ftype, ttype)

if strcmp(ttype, 'Diffusion')
    if strcmp(ftype, 'dfem')
        fcall = @perform_dfem_diffusion;
    elseif strcmp(ftype, 'cfem')
        fcall = @perform_cfem_diffusion;
    end
elseif strcmp(ttype, 'Transport')
    if strcmp(ftype, 'dfem')
        fcall = @perform_dfem_transport;
    elseif strcmp(ftype, 'cfem')
        fcall = @perform_cfem_transport;
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = solution_allocation(ndat, DoF)
ndof = DoF.TotalDoFs;
% Generate DoF space for each energy group and flux moment
% --------------------------------------------------------
if strcmp(ndat.transportMethod, 'Diffusion')
    mf = ndat.numberEnergyGroups;
    nf = 1;
else
    mf = ndat.numberEnergyGroups;
    nf = ndat.fluxMoments;
end
sol.flux = cell(mf,nf);
for i=1:mf
    for j=1:nf
        sol.flux{i,j} = ones(ndof,1);
    end
end
sol.numberEnergyGroups = mf;
sol.fluxMoments = nf;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = compute_flux_moment_differences(flux,flux0,ngs,mom)
global glob
err = 0;
for g=1:length(ngs)
    for m=1:length(mom)
        denom = sqrt(flux{g,m}'*flux{g,m});
        if denom < glob.small
            err = err + sqrt((flux{g,m} - flux0{g,m})'*(flux{g,m} - flux0{g,m}));
        else
            err = err + sqrt((flux{g,m} - flux0{g,m})'*(flux{g,m} - flux0{g,m}))/denom;
        end
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [keff,flux] = estimate_new_keff(flux,flux0)
ng = size(flux,1);
num = 0; denom = 0;
for g=1:ng
    num = num + norm(flux{g,1});
    denom = denom + norm(flux0{g,1});
end
keff = num/denom;
for g=1:ng
    flux{g,1} = flux{g,1} / keff;
end
return