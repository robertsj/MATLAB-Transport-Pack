%> @file  optimize_flare.m
%> @brief Demonstration optimization of FLARE parameters

clear classes

% Get the sample materials for fresh, once-, and twice-burned bundles.
mat = flare_mat();

% The bundle placement is given below.  Only quarter cores with rotational
% and albedo boundary conditions are supported.  Zeros must be places where
% a bundle is defined via rotational symmetry or where reflector exists.
% The bundle id's correspond to materials in the material database.
map = [   3  1  2  1  2  1  2  3  0
          0  2  1  2  1  2  1  3  0
          0  1  2  1  2  1  2  3  0
          0  2  1  2  1  2  3  3  0
          0  1  2  1  2  1  3  0  0
          0  2  1  2  1  2  3  0  0
          0  1  2  3  3  3  0  0  0
          0  3  3  3  0  0  0  0  0
          0  0  0  0  0  0  0  0  0];

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

refp=[  1.82459  2.38239  2.05785  2.09595  1.50127  1.22793  0.61444  0.1925   0. 
        2.38239  2.10583  2.2998   1.79482  1.65369  1.03418  0.66646  0.18749  0.    
        2.05785  2.2998   1.89146  1.88036  1.31323  1.02131  0.46492  0.1395   0.    
        2.09595  1.79482  1.88036  1.41645  1.24785  0.69612  0.26969  0.08316  0.    
        1.50127  1.65369  1.31323  1.24785  0.82396  0.57262  0.17873  0.       0.
        1.22793  1.03418  1.02131  0.69612  0.57262  0.28999  0.09297  0.       0.  
        0.61444  0.66646  0.46492  0.26969  0.17873  0.09297  0.       0.       0.   
        0.1925   0.18749  0.1395   0.08316  0.       0.       0.       0.       0.  
        0.       0.       0.       0.       0.       0.       0.       0.       0.]
    
