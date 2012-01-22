%> @file  test_fixed.m
%> @brief Tests the fixed source iteration
%
%> Finish me.
% ==============================================================================

clear classes

% Get the default input.
input = Input();

% Material etc.
put(input, 'number_groups',         2);

% Inner iteration parameters.
put(input, 'inner_tolerance',       1e-3);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  6);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');
% One material, one group, c = 0.9
mat         = test_materials(2);

% Simple uniformly meshed square
mesh        = test_mesh(2);

state       = State(input, mesh);
quadrature  = LevelSymmetric(4);
boundary    = Boundary(input, mesh, quadrature);

% Uniform source.
q_e         = Source(mesh, 2);
spectra     = [ 1.0 1
                0.0 0];      
placement   = [ 1 2; 2 2];    
set_sources(q_e, spectra, placement);

% Empty fission source.
q_f = FissionSource(state, mesh, mat);

% Make the inner iteration.
solver= Fixed(input,        ...
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
