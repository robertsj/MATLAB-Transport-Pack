%> @file  test_flare.m
%> @brief 2-d flare model.

%clear classes

% Get the sample materials for fresh, once-, and twice-burned bundles.
mat = flare_mat();

% The bundle placement is given below.  Only quarter cores with rotational
% and albedo boundary conditions are supported.  Zeros must be places where
% a bundle is defined via rotational symmetry or where reflector exists.
% The bundle id's correspond to materials in the material database.
% map = [   3  1  2  1  2  1  2  3  0
%           0  2  1  2  1  2  1  3  0
%           0  1  2  1  2  1  2  3  0
%           0  2  1  2  1  2  3  3  0
%           0  1  2  1  2  1  3  0  0
%           0  2  1  2  1  2  3  0  0
%           0  1  2  3  3  3  0  0  0
%           0  3  3  3  0  0  0  0  0
%           0  0  0  0  0  0  0  0  0];
map = [   1 1 1 0 
          0 1 1 0
          0 1 1 0
          0 0 0 0];
% Generate the input and add FLARE-specific things.  
input = Input();
put(input, 'node_width',    21);
put(input, 'mixing_factor', 0.88);
put(input, 'albedo_single', 0.22);
put(input, 'albedo_double', 0.65); 

% Create the solver.  This generates the coupling coefficients.
solver = FLARE(input, mat, map);

% Solve, returning the fission density and k-eigenvalue.
[s, k]=solve(solver);

% Plot the peaking factor.
plot_peak(solver);
