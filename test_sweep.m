function test_sweep()

tic

% Test input.
input = Input();

% Get mesh and materials
mesh = test_mesh(2);
mat  = test_materials(2);

% Quadrature
quadrature = LevelSymmetric(16);

% Boundary
boundary = Boundary(input, mesh, quadrature);
initialize(boundary, 1);

% Make DD equation
equation = DD2D(mesh, mat);
setup_group(equation, 1);

% Make sweeper
sweeper = Sweep2D(input, mesh, mat, quadrature, boundary, equation);
% Make a sweep source
source  = ones(number_cells(mesh), 1);
phi     = sweep(sweeper, source, 1);
plot_map(phi)

toc
end