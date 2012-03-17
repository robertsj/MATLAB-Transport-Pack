%> @file  test_eigenvalue.m
%> @brief Tests the eigenvalue solver in 2-D.
%
%> The tests are the same as the 1-D tests, using reflective conditions
%> to simulate 1-D.
% ==============================================================================

clear classes


% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-5);
put(input, 'eigen_max_iters',       100);
put(input, 'inner_tolerance',       0.4);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);


% ==============================================================================
% MESH 
% ==============================================================================
% Assembly coarse mesh edges
base   = [ 1.1580 4.4790 7.8000 11.1210 14.4420 15.6000 ]; 
% and fine mesh counts.
basef  = [ 1 2 2 2 2 1 ]*2; 
% Several such assemblies to make the total coarse mesh definition
xcm    = [ 0.0  base  base+15.6  base+15.6*2 base+15.6*3 ...
           base+15.6*4 base+15.6*5 base+15.6*6 ];
xfm    = [ basef basef basef basef basef basef basef ];

% Assembly types
assem_A    = [ 1 2 3 3 2 1 ];
assem_B    = [ 1 2 2 2 2 1 ];
assem_C    = [ 1 2 4 4 2 1 ];
assem_D    = [ 1 4 4 4 4 1 ];

% Cores 1, 2 and 3
core_1 = [ assem_A assem_B assem_A assem_B assem_A assem_B assem_A ];
core_2 = [ assem_A assem_C assem_A assem_C assem_A assem_C assem_A ];
core_3 = [ assem_A assem_D assem_A assem_D assem_A assem_D assem_A ];
ycm = [0 1];
yfm =  1;
mesh = Mesh2D(xfm, yfm, xcm, ycm, core_1);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = LevelSymmetric(16);
boundary    = BoundaryMesh(input, mesh, quadrature);
q_e         = Source(mesh, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, mesh, mat);  % Inititalized = used.
initialize(q_f);

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Eigensolver(input,         ...
                     state,         ...
                     boundary,      ...
                     mesh,          ...
                     mat,           ...
                     quadrature,	...
                     q_e,           ...
                     q_f);
 
tic
out = solve(solver); 
toc

% ==============================================================================
% POSTPROCESS 
% ==============================================================================

subplot(2, 1, 1)
f = flux(state, 1);
plot_flux(mesh, f), axis equal
subplot(2, 1, 2)
f = flux(state, 2);
plot_flux(mesh, f), axis equal
