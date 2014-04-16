function [data, geometry] = load_user_input()
global glob
% Neutronics Data
% ---------------
data.Neutronics.transportMethod = 'Diffusion';
data.Neutronics.FEMType = 'CFEM';
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

%%%%%%%%%%%%%%%%%%%%%%
% Transport Properties
%%%%%%%%%%%%%%%%%%%%%%
data.Neutronics.Transport.QuadType = 'LDFE';
data.Neutronics.Transport.angularQuadrature = 2;
data.Neutronics.Transport.fluxMoments = 1;
data.Neutronics.Transport.performDSA = 0;
% Physical Properties
data.Neutronics.Transport.TotalXS = [1.0];
data.Neutronics.Transport.AbsorbXS = [0.1];
data.Neutronics.Transport.ScatteringXS = [0.9];
data.Neutronics.Transport.FissionXS = [0.3];
data.Neutronics.Transport.FissSpec = [1.0];
data.Neutronics.Transport.ExtSource = [1.0];
% Boundary Conditions


%%%%%%%%%%%%%%%%%%%%%%
% Diffusion Properties
%%%%%%%%%%%%%%%%%%%%%%
% Physical Properties
data.Neutronics.Diffusion.DiffXS = [1/3];
data.Neutronics.Diffusion.TotalXS = [1.0];
data.Neutronics.Diffusion.AbsorbXS = [0.5];
data.Neutronics.Diffusion.ScatteringXS = [0.5];
data.Neutronics.Diffusion.FissionXS = [0.5];
data.Neutronics.Diffusion.FissSpec = [1.0];
data.Neutronics.Diffusion.ExtSource = [0.0];
% Boundary Conditions
data.Neutronics.Diffusion.BCFlags = [glob.Dirichlet];
data.Neutronics.Diffusion.BCVals = [0];

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-12;
data.solver.relativeTolerance = 1e-10;
data.solver.maxIterations = 2000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Problem Input Parameters
% ------------------------
% data.problem.problemType = 'SourceDriven';
data.problem.problemType = 'Eigenvalue';
data.problem.refineMesh = 0;
data.problem.refinementLevels = 0;
data.problem.refinementTolerance = 0.5;
data.problem.refinementType = 1;
data.problem.refinementSplits = 1;
data.problem.plotSolution = 0;

% Geometry Data
% -------------
dim = 2;
geometry = GeneralGeometry(dim, 'triangle', 'geometry_inputs/square.3');

