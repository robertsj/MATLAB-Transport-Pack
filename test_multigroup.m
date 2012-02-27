%> @file  test_multigroup.m
%> @brief Tests and compares Krylov to Gauss-Seidell for multigroup.
% ==============================================================================

%clear classes

input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_tolerance',       1e-8);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-8);
put(input, 'outer_max_iters',       300);
put(input, 'inner_solver',          'GMRES');
put(input, 'bc_left',               'vacuum');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'vacuum');
put(input, 'bc_bottom',             'vacuum');
put(input, 'print_out',             1);

% Two materials, two groups, some upscatter
mat         = test_materials(3);

% Simple uniformly meshed square
mesh        = test_mesh(3);

state       = State(input, mesh);
quadrature  = QuadrupleRange(8);

boundary    = BoundaryMesh(input, mesh, quadrature);

% Empty external source.
q_e         = Source(mesh, 2);
spectra     = [ 1.0 1.0
                1.0 1.0];      
placement   = [ 2 2; 2 2];    
set_sources(q_e, spectra, placement);

% Use a fission source.
q_f = FissionSource(state, mesh, mat);
%initialize(q_f);

% Make the inner iteration.
solver= KrylovMG(input,        ...
              state,    	...
              boundary,     ...
              mesh,     	...
              mat,        	...
              quadrature, 	...
              q_e,          ...
              q_f);
 % Solve the problem
tic
out = solve(solver); 
toc

% Get the flux
f = flux(state, 1);
subplot(2,1,1)
plot_flux(mesh, f)
axis square
% Get the flux
f = flux(state, 2);
% and plot it
subplot(2,1,2)
plot_flux(mesh, f)
axis square
