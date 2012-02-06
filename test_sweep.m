

tic

% Test input.
input = Input();

% Get mesh and materials
mesh = test_mesh(2);
mat  = test_materials(2);

% Quadrature
quadrature = LevelSymmetric(16);

% Boundary
boundary = BoundaryMesh(input, mesh, quadrature);
initialize(boundary, 1);

% Make DD equation
equation = DD2D(mesh, mat, quadrature);
setup_group(equation, 1);

% Make sweeper
sweeper = Sweep2D_mod(input, mesh, mat, quadrature, boundary, equation);
% Make a sweep source
source  = ones(number_cells(mesh), 1);
phi     = sweep(sweeper, source, 1);
plot_flux(mesh, phi)

toc
