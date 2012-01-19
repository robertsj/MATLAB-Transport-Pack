%> @file  test_eigenvalue.m
%> @brief Tests the eigenvalue solver.
%
%> Finish me.
% ==============================================================================

clear


% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_solver',          'Livolant');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);


% ==============================================================================
% MESH 
% ==============================================================================
xcm    = [ 0.0  45.0  50.0];
xfm    = [    45     5    ];
ycm    = [ 0.0  45.0  50.0];
yfm    = [    45     5    ];
mat_map = [ 1  1     % ^
            2  1 ];  % | y   x -->
mesh = Mesh2D(xfm, yfm, xcm, ycm, mat_map);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = LevelSymmetric(2);
boundary    = Boundary(input, mesh, quadrature);
q_e         = Source(mesh, 2);          % Not initialized = not used.
q_f = FissionSource(state, mesh, mat);  % Inititalized = used.
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
plot_flux(mesh, f)
subplot(2, 1, 2)
f = flux(state, 2);
plot_flux(mesh, f)
