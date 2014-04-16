function out = setGlobals ()

% Numbers
%%%%%%%%%
out.small = 1e-14;
out.large = 1/eps/100;

% Constants
%%%%%%%%%%%
out.speed_of_light = 2.998e10;
out.planck = 6.626176e-34;
out.boltzmann_ev = 8.617e-5;
out.boltzmann_J = 1.380662e-23;
out.elementary_charge = 1.6021892e-19;
out.neutron_mass = 1.6749544e-27;
out.avogadro = 6.022e23;
out.gas_constant = 8.31441;

out.amu_to_kg = 1.6605655e-27;
out.amu_to_MeV = 931.5016;
out.MeV_to_J = 1.601892e-13;
out.J_to_MeV = 6.2415e12;

% Miscellaneous
%%%%%%%%%%%%%%%
out.input_path = '..\Input\';
out.output_path = '..\Output\';
out.directory = [];

% Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%
out.Dirichlet = 1;
out.Neumann = 2;
out.Robin = 3;

% Transport Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.Vacuum = 1;
out.Reclecting = 2;
out.IncidentIsotropic = 3;
out.IncidentCurrent = 4;