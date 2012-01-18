%> @file  test_equation.m
%> @brief Tests the equation classes
function test_equation

% Get mesh and materials
mesh = test_mesh();
mat  = test_materials();

% Make DD equation
equation = DD2D(mesh, mat);

% Setup.
setup_angle(equation, 0.33, 0.66);

% Solve
psi_in = [0.25; 0.35];
solve(equation, psi_in, 1.0, 1, 1, 1)

end